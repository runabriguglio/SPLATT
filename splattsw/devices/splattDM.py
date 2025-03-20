import numpy as np

class SPLATTEngine():

    def __init__(self, ip:str = '193.206.155.220', port:int = 9090):
        import Pyro4
        self.eng = Pyro4.Proxy(f"PYRO:matlab_engine@{ip}:{port}")
        self.eng.connect_matlab()

        print('Initializing mirror variables...')
        self._shellset = False
        self.nActs = int(self.eng.read_data('sys_data.mirrNAct'))
        self.actCoords = self._get_act_coords()
        self._bits2meters = float(self.eng.read_data('2^-sys_data.coeffs.Scale_F_Lin'))


    def get_position(self):
        if self.shellset is False:
            print('Shell must be set before giving commands!')
            self._set_shell()
        posCmdBits = np.array(self.eng.read_data("aoRead('sabu16_position',1:19)"))
        posCmd = posCmdBits * self._bits2meters
        return posCmd


    def set_position(self, cmd):
        if self.shellset is False:
            print('Shell must be set before giving commands!')
            self._set_shell()
        self.eng.send_command(f"splattMirrorCommand({cmd},'relative')")


    def read_buffers(self, n_samples:int = 128, decimation:int = 0):

        if n_samples > 256:
            raise ValueError('Maximum number of samples is 256!')

        self.send_command(f"clear opts; opts.dec = {decimation}; opts.sampleNr = {n_samples}; opts.save2fits = 1; opts.save2mat = 0; opts.saveCmds = 1")
        self.send_command("[pos,cur,buf_tn]=splattAcqBufInt({'sabi32_Distance','sabi32_pidCoilOut'},opts)")

        buf_tn = self.eng.read_data('buf_tn')
        mean_pos = np.array(self.eng.read_data('mean(pos,2)')*self._bits2meters)
        mean_cur = np.array(self.eng.read_data('mean(cur,2)'))

        return mean_pos, mean_cur, buf_tn


    def _set_shell(self):

        try:
            pos = np.array(self.eng.read_data('lattGetPos()'))
            rest_pos = np.array(self.eng.read_data('sys_data.restpos'))
            if min(pos) > max(rest_pos):
                self._shellset = True
        except:
            print('Performing startup ...')
            self.eng.send_command('splattStartup')
            
        if self._shellset is False:
            print('Setting the shell...')
            self.eng.send_command('splattFastSet')
            self._shellset = True

    def _get_act_coords(self):
        self.eng.send_command("phi = deg2rad(60)")
        self.eng.send_command("matrixRot = [cos(phi),-sin(phi);sin(phi),cos(phi)]")
        self.eng.send_command("coord_Act = matrixRot*mirrorData.coordAct'")
        coordAct = np.array(self.eng.read_data("coord_Act'"))
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
        shape = self._dm.get_position()
        return shape

    def set_shape(self, cmd, differential:bool=False):
        if differential:
            shape = self._dm.get_position()
            cmd = cmd + shape
        self._checkCmdIntegrity(cmd)
        self._dm.set_position(cmd)

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

