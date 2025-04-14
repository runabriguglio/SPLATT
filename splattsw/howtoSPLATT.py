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
from m4.dmutils.flattening import Flattening

'''
icp = IFFCapturePreparation(dm)
cmh = icp.createTimedCmdHistory(modesList=[0,1,2,3,4,5,6],modesAmp=1e-6)
dm.uploadCmdHistory(cmh)
dm.runCmdHistory(interf)
'''
#the following command makes all together
mlist = [0,1,2,3,4,5,6]
mamp = 5e-6
nmodes2flat = len(mlist)
nmodes2discard = 3
tn = ifm.iffDataAcquisition(dm, interf,mlist,amplitude=mamp)

rebinfact = 4
#adm.set_shape(np.zeros(19))
ifp.process(tn,  save_and_rebin_cube=(True,rebinfact))

fl = Flattening(tn)
fl.filterIntCube([1,2,3])
fl.applyFlatCommand(dm, interf, mlist, nframes=1,modes2discard=nmodes2discard)

img = interf.acquire_map(1, rebin=rebinfact)
#fl.applyFlatCommand(adm, interf, 820, modes2discard=modes2remove[x])
fl.loadImage2Shape(img)
fl.computeRecMat(nmodes2discard)
deltacmd = fl.computeFlatCmd(nmodes2flat)
dm.set_shape(deltacmd, differential=True)
time.sleep(1)
fimg = interf.acquire_map(5, rebin=4))

from splattsw import acceleration_analysis as acc
from splattsw.devices.webDAQ import WebDAQ

wbdq = WebDAQ()
wbdq.connect()

ffv = dm.mirrorModes
mode_id = 4
mode_amp = 5e-6
cmd = ffv[:,mode_id]*mode_amp

wbdq.start_schedule()
tn_buf = dm.sendBuffer))
tn = interf.capture(500)
wbdq.stop_schedule()

acc.wdsync()
wdfile = acc.last_wdfile()
interf.produce(tn)


from splattsw import userscripts as usr
import matplotlib.pyplot as plt

v = usr.analyze_opt_step(tn)
em = usr.analyze_buf_step(tn_buf, ffv)
tvec = em['time']
pos = em['position']
pos_cmd = em['position_command']
plt.figure()
plt.plot(tvec,pos[:,1:])
plt.plot(tvec,pos_cmd[:,1:],'--')
plt.grid(True)
plt.title(tn_buf)

pos_rms = em['pos_rms']
plt.figure()
plt.plot(tvec,pos_rms,label='position RMS (first 4 modes removed)')
plt.grid(True)
plt.legend()
plt.title(tn_buf)

interf_freq = 80
tt = np.arange(len(v))/interf_freq
plt.figure()
plt.plot(tt,v)
plt.grid(True)
plt.title(tn)

data = acc.openfile(wdfile)
freq_WebDAQ = 1651
tacc = np.arange(np.shape(data)[-1])/freq_WebDAQ
plt.figure()
plt.plot(tacc,data[0:2].T)
plt.legend('Channel 0','Channel 1')
plt.title(wdfile)

