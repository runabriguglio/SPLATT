
class SPLATTEngine():

    def __init__(self, ip:str = '193.206.155.220', port:int = 9090):
        import Pyro4
        self.eng = Pyro4.Proxy(f"PYRO:matlab_engine@{ip}:{port}")
        print('Starting the matlab engine...')
        self.eng.start_engine()

        print('Performing startup...')
        self.eng.send_command('splattInit')
        #self.eng.send_command('splattStartup')

        print('Initializing mirror variables...')
        self._shellIsSet = False
        self.nActs = int(self.eng.get_data('sys_data.mirrNAct'))
        self.actCoords = self._get_act_coords()
        self._bits2meters = float(self.eng.get_data('2^-sys_data.coeffs.Scale_F_Lin'))


    def get_shape(self):
        if self._shellIsSet is False:
            print('Shell must be set before giving commands!')
            self._set_shell()
        posCmdBits = np.array(self.eng.get_data("aoRead('sabu16_position',1:19)"))
        posCmd = posCmdBits * self._bits2meters
        return posCmd


    def set_shape(self, cmd):
        if self._shellIsSet is False:
            print('Shell must be set before giving commands!')
            self._set_shell()
        self.eng.send_command(f"splattMirrorCommand(cmd,'relative')")


    def _set_shell(self):
        print('Setting the shell...')
        self.eng.send_command('splattFastSet')
        self._shellset = True

    def _get_act_coords(self):
        self.eng.send_command("phi = deg2rad(60)")
        self.eng.send_command("matrixRot = [cos(phi),-sin(phi);sin(phi),cos(phi)]")
        self.eng.send_command("coord_Act = matrixRot*mirrorData.coordAct'")
        coordAct = np.array(self.eng.get_data("coord_Act'"))
        return coordAct



class SPLATTDm(BaseDeformableMirror):
    """
    SPLATT interface with M4 software.
    """

    def __init__(self, ip:str = '193.206.155.220', port:int = 9090):
        """The Constructor"""
        self._dm            = SPLATTEngine(ip,port)
        self.nActs          = self._dm.nActs
        self.mirrorModes    = None
        self.actCoord       = self._dm.actCoords
        self.cmdHistory     = None
        self.baseDataPath   = fn.OPD_IMAGES_ROOT_FOLDER
        self.refAct         = 1

    def get_shape(self):
        shape = self._dm.get_shape()
        return shape

    def set_shape(self, cmd, differential:bool=False):
        if differential:
            shape = self._dm.get_shape()
            cmd = cmd + shape
        self._checkCmdIntegrity(cmd)
        self._dm.set_shape(cmd)

    def uploadCmdHistory(self, cmdhist):
        self.cmdHistory = cmdhist

    def runCmdHistory(self, interf=None, delay=2., save:str=None, differential:bool=True):
        if self.cmdHistory is None:
            raise ValueError("No Command History to run!")
        else:
            tn = _ts.now() if save is None else save
            print(f"{tn} - {self.cmdHistory.shape[-1]} images to go.")
            datafold = os.path.join(self.baseDataPath, tn)
            s = self.get_shape()
            if not os.path.exists(datafold) and interf is not None:
                os.mkdir(datafold)
            for i,cmd in enumerate(self.cmdHistory.T):
                print(f"{i+1}/{self.cmdHistory.shape[-1]}", end="\r", flush=True)
                if differential:
                    cmd = cmd+s
                self.set_shape(cmd)
                if interf is not None:
                    time.sleep(delay)
                    img = interf.acquire_phasemap()
                    path = os.path.join(datafold, f"image_{i:05d}.fits")
                    rd.save_phasemap(path, img)
        self.set_shape(s)
        return tn


    def nActuators(self):
        return self.nActs

    def _checkCmdIntegrity(self, cmd):
        mcmd = np.max(cmd)
        if mcmd > 2e-6:
            raise ValueError(f"Command value {mcmd} is greater than 2 [um]")
        mcmd = np.min(cmd)
        if mcmd < -2e-6:
            raise ValueError(f"Command value {mcmd} is smaller than - 2 [um]")
        scmd = np.std(cmd)
        if scmd > 0.1:
            raise ValueError(f"Command standard deviation {scmd} is greater than 0.1.")

