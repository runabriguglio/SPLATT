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

     def set_shape(self, cmd, add2LastCmd:bool=False):
         if add2LastCmd:
             cmd = cmd + self.get_shape()
         cmd = cmd - self._dm.flatPos # set_position is relative to flatPos
         self._checkCmdIntegrity(cmd) # now check integrity on the relative command
         self._dm.set_position(cmd) 

     def uploadCmdHistory(self, cmdhist):
         self.cmdHistory = cmdhist

     def runCmdHist(self, interf=None, delay=0.2, save_tn:str=None, differential:bool=True):
         if self.cmdHistory is None:
             raise ValueError("No Command History to run!")
         else:
             tn = _ts.now() if save_tn is None else save_tn
             print(f"{tn} - {self.cmdHistory.shape[-1]} images to go.")
             datafold = os.path.join(self.baseDataPath, tn)
             s = self._dm.flatPos # self.get_shape()
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

     def nActuators(self):
         return self.nActs
     
     def saveFlatTN(self):
         tn = self._dm._eng.read_data('tn=lattSaveFlat()')
         return tn
     
     def loadFlatTN(self, tn:str):
         self._dm._eng.send_command(f"lattLoadFlat('{tn}')")
         self._dm.flatPos, self._dm.flatForce = self._dm.read_flat_data()

     def _checkCmdIntegrity(self, cmd):
         pos = cmd + self._dm.flatPos
         if np.max(pos) > 1.2e-3:
            raise ValueError(f'End position is too high at {np.max(pos)*1e+3:1.2f} [mm]')
         if np.min(pos) < 450e-6:
            raise ValueError(f'End position is too low at {np.min(pos)*1e+3:1.2f} [mm]')
         # Check force
         force = self._dm.flatForce + self._dm.ffMatrix * (cmd - np.mean(cmd))
         force_thr = self._dm.maxForce
         if np.max(np.abs(force)) >= force_thr:
             raise ValueError(f"Command would require {force_thr} [N]")


class SPLATTEngine():

    def __init__(self, ip:str = '193.206.155.220', port:int = 9090):

        import Pyro4
        self._eng = Pyro4.Proxy(f"PYRO:matlab_engine@{ip}:{port}")
        
        print('Reading mirror variables ...')
        self.nActs = int(self._eng.read_data('sys_data.mirrNAct'))
        self.actCoords = np.array(self._eng.read_data('mirrorData.coordAct'))
        self.mirrorModes = np.array(self._eng.read_data('sys_data.ff_v'))
        self.ffMatrix = np.array(self._eng.read_data('sys_data.ff_matrix'))

        self._bits2meters = float(self._eng.read_data('2^-sys_data.coeffs.Scale_F_Lin'))
        self._N2bits = float(self._eng.read_data('sys_data.coeffs.Force2DAC_V'))
        self._satThreshold = float(self._eng.read_data('sys_data.currentSatThreshold'))

        self._shellset = True
        try:
            self.flatPos, self.flatForce = self.read_flat_data()
        except:
            self._shellset = False
            print('Unable to read set position: remember to perform startup and set shell')

        self.maxForce = self._satThreshold / self._N2bits

    def get_position(self):
        posCmdBits = np.array(self._eng.read_data("aoRead('sabu16_position',1:19)"))
        posCmd = posCmdBits * self._bits2meters
        posCmd = np.reshape(posCmd, self.nActs)
        return posCmd

    def set_position1(self, cmd): 
        if self._shellset is False: print('Shell must be set before giving commands!')
        cmd = cmd.tolist()
        self._eng.send_command(f"splattMirrorCommand({cmd}','zonal')")

    def set_position2(self, cmd):
        if self._shellset is False: print('Shell must be set before giving commands!')
        cmd += self.flatPos
        cmd_bits = cmd / self._bits2meters
        cmd_bits = cmd_bits.tolist()
        f_cmd = self.ffMatrix * cmd
        f_cmd_bits = f_cmd * self._N2bits + self.flatForce
        f_cmd_bits = f_cmd_bits.tolist()
        self._eng.send_command(f"aoWrite('sabu16_position',{cmd_bits}',1:19)")
        self._eng.send_command(f"lattApplyForce({f_cmd_bits}')")

    def set_position(self, cmd):
        if self._shellset is False: print('Shell must be set before giving commands!')
        cmd += self.flatPos
        cmd_bits = cmd / self._bits2meters
        cmd_bits = cmd_bits.tolist()
        self._eng.send_command(f"aoWrite('sabu16_position',{cmd_bits}',1:19)")


    def read_buffers(self, external: bool = False, n_samples:int = 128, decimation:int = 0):

        if n_samples > 256:
            raise ValueError('Maximum number of samples is 256!')

        self._eng.send_command(f"clear opts; opts.dec = {decimation}; opts.sampleNr = {n_samples}; opts.save2fits = 1; opts.save2mat = 0")
        print('Reading buffers, hold tight: this may take a while ...')
        if external:
            self._eng.send_command("[pos,cur,buf_tn]=splattAcqBufExt({'sabi32_Distance','sabi32_pidCoilOut'},opts)")
        else:
            self._eng.send_command("[pos,cur,buf_tn]=splattAcqBufInt({'sabi32_Distance','sabi32_pidCoilOut'},opts)")

        buf_tn = self._eng.read_data('buf_tn')
        mean_pos = np.array(self._eng.read_data('mean(pos,2)')) * self._bits2meters
        mean_pos = np.reshape(mean_pos,self.nActs)
        mean_cur = np.array(self._eng.read_data('mean(cur,2)'))
        mean_cur = np.reshape(mean_cur, self.nActs)

        return mean_pos, mean_cur, buf_tn
    

    def read_flat_data(self):
        flatPos = np.array(self._eng.read_data('sys_data.flatPos')) * self._bits2meters
        flatPos = np.reshape(flatPos,self.nActs)
        flatForce = np.array(self._eng.read_data('sys_data.flatCur')) / self._N2bits
        flatForce = np.reshape(flatForce,self.nActs)
        return flatPos, flatForce

    def _set_shell(self):

        if self._shellset is False:
            print('Setting the shell...')
            self._eng.send_command('splattFastSet')
            self._shellset = True
        else:
            print('Shell set variable is True, overwrite it if you wish to set again')

