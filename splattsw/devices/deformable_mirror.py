"""
Author(s)
    - Pietro Ferraiuolo : written in 2024
    - Runa Briguglio : written in 2024

Description
-----------
This module contains the class that defines the M4 Adaptive Unit (M4AU) device.
"""
import os
import time
import numpy as np
from astropy.io import fits as pyfits
from m4.configuration import config_folder_names as fn
from m4.devices.base_deformable_mirror import BaseDeformableMirror
from m4.ground import logger_set_up as lsu, timestamp, read_data as rd

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
             shape = self._dm.get_position()
             cmd = cmd + shape
         self._checkCmdIntegrity(cmd)
         self._dm.set_position(cmd)

     def uploadCmdHistory(self, cmdhist):
         self.cmdHistory = cmdhist

     def runCmdHistory(self, interf=None, delay=0.5, save:str=None, differential:bool=True):
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
         if mcmd > 5e-6:
             raise ValueError(f"Command value {mcmd} is greater than 5 [um]")
         mcmd = np.min(cmd)
         if mcmd < -5e-6:
             raise ValueError(f"Command value {mcmd} is smaller than -5 [um]")
         scmd = np.std(cmd)
         if scmd > 0.1:
             raise ValueError(f"Command standard deviation {scmd} is greater than 0.1.")


class SPLATTEngine():

    def __init__(self, ip:str = '193.206.155.220', port:int = 9090):

        import Pyro4
        self._eng = Pyro4.Proxy(f"PYRO:matlab_engine@{ip}:{port}")

        print('Initializing mirror variables...')
        self.nActs = int(self._eng.read_data('sys_data.mirrNAct'))
        self.actCoords = np.array(self._eng.read_data('mirrorData.coordAct'))
        self.mirrorModes = np.array(self._eng.read_data('sys_data.ff_v'))

        self._shellset = False
        self._bits2meters = float(self._eng.read_data('2^-sys_data.coeffs.Scale_F_Lin'))


    def get_position(self):
        if self.shellset is False:
            print('Shell must be set before giving commands!')
            self._set_shell()
        posCmdBits = np.array(self._eng.read_data("aoRead('sabu16_position',1:19)"))
        posCmd = posCmdBits * self._bits2meters
        return posCmd


    def set_position(self, cmd):
        if self.shellset is False:
            print('Shell must be set before giving commands!')
            self._set_shell()
        self._eng.send_command(f"splattMirrorCommand({cmd},'relative')")


    def get_ff_force(self):
        ff_force = np.array(self._eng.read_data("aoRead('sabi16_force',1:19)"))
        return ff_force
    

    def read_buffers(self, n_samples:int = 128, decimation:int = 0):

        if n_samples > 256:
            raise ValueError('Maximum number of samples is 256!')

        self.send_command(f"clear opts; opts.dec = {decimation}; opts.sampleNr = {n_samples}; opts.save2fits = 1; opts.save2mat = 0")
        print('Reading buffers, hold tight ...')
        self.send_command("[pos,cur,buf_tn]=splattAcqBufInt({'sabi32_Distance','sabi32_pidCoilOut'},opts)")

        buf_tn = self._eng.read_data('buf_tn')
        mean_pos = np.array(self._eng.read_data('mean(pos,2)')*self._bits2meters)
        mean_cur = np.array(self._eng.read_data('mean(cur,2)'))

        return mean_pos, mean_cur, buf_tn


    def _set_shell(self):

        try:
            pos = np.array(self._eng.read_data('lattGetPos()'))
            rest_pos = np.array(self._eng.read_data('sys_data.restpos'))
            if min(pos) > max(rest_pos):
                self._shellset = True
        except:
            print('Performing startup ...')
            self._eng.send_command('splattStartup')

        if self._shellset is False:
            print('Setting the shell...')
            self._eng.send_command('splattFastSet')
            self._shellset = True

