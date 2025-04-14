"""
Author(s)
    - Pietro Ferraiuolo : written in 2024
    - Runa Briguglio : written in 2024
    - Matteo Menessini : adapted for SPLATT in 2025

Description
-----------
This module contains the class that defines the SPLATT deformable mirror device.
"""
import os
import time
import numpy as np
from astropy.io import fits as pyfits
from m4.configuration import config_folder_names as fn
from m4.devices.base_deformable_mirror import BaseDeformableMirror
from m4.ground import logger_set_up as lsu, timestamp, read_data as rd

_ts = timestamp.Timestamp()

class SPLATTDm(BaseDeformableMirror):
     """
     SPLATT interface with M4 software.
     """

     def __init__(self, ip:str = '193.206.155.220', port:int = 9090):
         """The Constructor"""
         self._dm            = SPLATTEngine(ip,port)
         self.nActs          = self._dm.nActs
         self.mirrorModes    = self._dm.mirrorModes
         self.actCoord       = self._dm.actCoords
         self.cmdHistory     = None
         self.baseDataPath   = fn.OPD_IMAGES_ROOT_FOLDER
         self.refAct         = 16

     def get_shape(self):
        shape = self._dm.get_position()
        return shape

     def set_shape(self, cmd, differential:bool=False):
         if differential:
            lastCmd = self._dm.get_position_command()
            cmd = cmd + lastCmd
         self._checkCmdIntegrity(cmd)
         self._dm.set_position(cmd) 

     def uploadCmdHistory(self, cmdhist):
         self.cmdHistory = cmdhist

     def runCmdHistory(self, interf=None, delay=0.2, save:str=None, differential:bool=True):
         if self.cmdHistory is None:
             raise ValueError("No Command History to run!")
         else:
             tn = _ts.now() if save is None else save
             print(f"{tn} - {self.cmdHistory.shape[-1]} images to go.")
             datafold = os.path.join(self.baseDataPath, tn)
             s = self._dm.get_position_command()  #self._dm.flatPos # self.get_shape()
             if not os.path.exists(datafold) and interf is not None:
                 os.mkdir(datafold)
             for i,cmd in enumerate(self.cmdHistory.T):
                 print(f"{i+1}/{self.cmdHistory.shape[-1]}", end="\r", flush=True)
                 if differential:
                     cmd = cmd+s
                 self.set_shape(cmd)
                 if interf is not None:
                     time.sleep(delay)
                     img = interf.acquire_map()
                     path = os.path.join(datafold, f"image_{i:05d}.fits")
                     rd.save_phasemap(path, img)
         self.set_shape(s)
         return tn

    def get_forces(self):
        forces = self._dm.get_force()
        return forces

     def sendBufferCommand(self, cmd, differential:bool=False, cmdPreTime:float = 10e-3, freq:float = None, delay = 1.0): 
        # cmd is a command relative to self._dm.flatPos
         if differential:
            lastCmd = self._dm.get_position_command()
            cmd = cmd + lastCmd
         self._checkCmdIntegrity(cmd) 
         cmd = cmd.tolist()
         if freq is not None:
            tn = self._dm._eng.read(f'prepareDynCmdHistory({cmd},{freq})')   
         else:
            tn = self._dm._eng.read(f'prepareCmdHistory({cmd},{cmdPreTime})')
         self._dm._eng.oneway_send(f'pause({delay}); sendCmdHistory(buffer)')
         return tn

     def nActuators(self):
         return self.nActs

     def _checkCmdIntegrity(self, cmd):
         pos = cmd + self._dm.flatPos
         if np.max(pos) > 1.2e-3:
            raise ValueError(f'End position is too high at {np.max(pos)*1e+3:1.2f} [mm]')            tn = self._dm._eng.read(f'prepareCmdHistory({cmd},{cmdPreTime})')

         if np.min(pos) < 450e-6:
            raise ValueError(f'End position is too low at {np.min(pos)*1e+3:1.2f} [mm]')


class SPLATTEngine():

    def __init__(self, ip:str = '193.206.155.220', port:int = 9090):

        import Pyro4
        self._eng = Pyro4.Proxy(f"PYRO:matlab_engine@{ip}:{port}")
        
        self.nActs = int(self._eng.read('sys_data.mirrNAct'))
        self.actCoords = np.array(self._eng.read('mirrorData.coordAct'))
        self.mirrorModes = np.array(self._eng.read('sys_data.ff_v'))
        self.ffMatrix = np.array(self._eng.read('sys_data.ff_matrix'))

        self._bits2meters = float(self._eng.read('2^-sys_data.coeffs.Scale_F_Lin'))
        self._N2bits = float(self._eng.read('sys_data.coeffs.Force2DAC_V'))
        
        self._shellset = True
        try:
            self.flatPos = self.read_flat_data()
        except:
            self._shellset = False
            print('Unable to read set position: remember to perform startup and set shell')

        print('Initialized SPLATT deformable mirror')

    def get_position_command(self): # relative to flatPos
        posCmdBits = np.array(self._eng.read("aoRead('sabu16_position',1:19)"))
        posCmd = posCmdBits * self._bits2meters
        posCmd = np.reshape(posCmd, self.nActs)
        posCmd -= self.flatPos
        return posCmd

    def get_position(self):
        pos = np.array(self._eng.read("lattGetPos()"))
        pos = np.reshape(pos, self.nActs)
        return pos

    def get_force(self):
        force = np.array(self._eng.read("lattGetForce()"))
        force = np.reshape(force, self.nActs)
        return force

    def set_position(self, cmd): 
        if self._shellset is False: print('Shell must be set before giving commands!')
        cmd = cmd.tolist()
        self._eng.send(f"splattMirrorCommand({cmd}')")


    def read_buffers(self, external: bool = False, n_samples:int = 128, decimation:int = 0):

        if n_samples > 256:
            raise ValueError('Maximum number of samples is 256!')

        self._eng.send(f"clear opts; opts.dec = {decimation}; opts.sampleNr = {n_samples}; opts.save2fits = 1; opts.save2mat = 0")
        print('Reading buffers, hold tight: this may take a while ...')
        if external:
            self._eng.send("[pos,cur,buf_tn]=splattAcqBufExt({'sabi32_Distance','sabi32_pidCoilOut'},opts)")
        else:
            self._eng.send("[pos,cur,buf_tn]=splattAcqBufInt({'sabi32_Distance','sabi32_pidCoilOut'},opts)")

        buf_tn = self._eng.read('buf_tn')
        mean_pos = np.array(self._eng.read('mean(pos,2)')) * self._bits2meters
        mean_pos = np.reshape(mean_pos,self.nActs)
        mean_cur = np.array(self._eng.read('mean(cur,2)'))
        mean_cur = np.reshape(mean_cur, self.nActs)

        return mean_pos, mean_cur, buf_tn
    
    def saveFlatTN(self, tn:str = None):
         if tn is None:
            tn = self._eng.read('lattSaveFlat()')
         else:
            tn = self._eng.read(f'lattSaveFlat({tn})')
         return tn

    def updateFlatTN(self, tn:str = None):
         if tn is not None:
            self._eng.send(f"lattLoadFlat('{tn}')")
         self.flatPos = self.read_flat_data()

    def read_flat_data(self):
        flatPos = np.array(self._eng.read('sys_data.flatPos')) * self._bits2meters
        flatPos = np.reshape(flatPos,self.nActs)
        return flatPos

    def _set_shell(self):

        if self._shellset is False:
            print('Setting the shell...')
            self._eng.send('splattFastSet')
            self._shellset = True
        else:
            print('Shell set variable is True, overwrite it if you wish to set again')

