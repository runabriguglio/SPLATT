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
from m4.dmutils import iff_module as ifm
from m4.dmutils.flattening import Flattening


icp = IFFCapturePreparation(dm)
cmh = icp.createTimedCmdHistory(modesList=[0,1,2,3,4,5,6],modesAmp=1e-6)
dm.uploadCmdHistory(cmh)
dm.runCmdHistory(interf)
#the following command makes all together
mlist = [0,1,2,3]
nmodes2flat = len(mlist)
nmodes2discard = 3
tn = ifm.iffDataAcquisition(dm, interf,mlist,amplitude=1e-6)

rebinfact = 4
#adm.set_shape(np.zeros(19))
ifp.process(tn,  save_and_rebin_cube=(True,rebinfact))

fl = Flattening(tn)
fl.filterIntCube([1,2,3])
img = interf.acquire_map(1, rebin=rebinfact)
#fl.applyFlatCommand(adm, interf, 820, modes2discard=modes2remove[x])
fl.loadImage2Shape(img)
fl.computeRecMat(nmodes2discard)
deltacmd = fl.computeFlatCmd(nmodes2flat)
dm.set_shape(deltacmd, differential=True)
time.sleep(1)
fimg = interf.acquire_phasemap(5, rebin=4))

ffv = dm.mirrorModes
mode_id = 7
mode_amp = 4e-6
cmd = ffv[:,mode_if]*mode_amp
tn_buf = dm.sendBufferCommand(cmd)
tn = interf.capture(500)


from splattsw import userscripts as usr

v = usr.analyze_opt_step(tn)




'''
