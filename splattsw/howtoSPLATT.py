#modules already imported by initSPLATT
'''
from m4 import main, noise
import os, time, numpy as np
from m4.devices.i4d import I4D
from m4.utils import osutils as osu
#from m4.mini_OTT import measurements
from matplotlib import pyplot as plt
#from m4.userscripts import OTTScripts
from m4.analyzers import timehistory as th
from m4.configuration.start import create_ott
from m4.ground import read_data as rd, zernike as zern
from m4.configuration import update_folder_paths as ufp
from splattsw.devices.deformable_mirror import SPLATTDm
import gentle
'''
#devices available (after initSPLATT)
'''
interf = gentle.PhaseCam()
dm = SPLATTDm()
'''

#IFF and flattening
from m4.dmutils import iff_module as ifm
from m4.dmutils.iff_acquisition_preparation import IFFCapturePreparation
from m4.dmutils import iff_processing as ifp


icp = IFFCapturePreparation(dm)
cmh = icp.createTimedCmdHistory(modesList=[0,1,2,3,4,5,6],modesAmp=1e-6)
dm.uploadCmdHistory(cmh)
dm.runCmdHist(interf)
#adm.set_shape(np.zeros(19))
ifp.process(tn, rebin=4, save_cube=True)




'''
