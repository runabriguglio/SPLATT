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
#general imports
##### launch      calpy -f /mnt/libero/SPLATTData/Data/SysConfig     #####
#### calpy is the bash command ---- equivalent to pyott  --- to create the environment variable and does the imports

#import aoptics
#pyconf = '/mnt/libero/SPLATTData/Data/SysConfigurations/configuration.yaml'
#aoptics.load_configuration_file(pyconf)

from opticalib.devices.interferometer import PhaseCam
from opticalib.devices.deformable_mirrors import SplattDm
from splattsw import splatt_analysis as sp
from opticalib.dmutils import iff_module as ifm
from opticalib.dmutils import iff_processing as ifp
from opticalib.dmutils.flattening import Flattening
import numpy as np
interf = PhaseCam('4020')

from splattsw.devices import moxa_io as mx
moxa=mx.Moxa_pt0()
from splattsw import measurements as me
meas = me.Measurements(interf)
tn = meas.opticalMonitoring(nframes, delay)

interf = PhaseCam('4020')
dm = SplattDm()

#mamp = 2e-6
mlist = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
mamp = np.array([2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1.5,1.5,1.5])*1e-6

nmodes2flat = len(mlist)
nmodes2discard = 3
tn = ifm.iffDataAcquisition(dm, interf,mlist,amplitude=mamp)

rebinfact = 4
ifp.process(tn,  save=True, rebin=rebinfact)

#tn = '20251008_152616'
fl = Flattening(tn)
fl.filterIntCube([1,2,3])
#fl.applyFlatCommand(dm, interf, mlist, nframes=1,modes2discard=nmodes2discard)

nimgs = 1
img = interf.acquire_map(nimgs, rebin=rebinfact)
#fl.applyFlatCommand(adm, interf, 820, modes2discard=modes2remove[x])
fl.loadImage2Shape(img)
fl.computeRecMat(nmodes2discard)
deltacmd = fl.computeFlatCmd(nmodes2flat)
cmdOffset = deltacmd
dm.set_shape(deltacmd, differential=True)
time.sleep(1)
fimg = interf.acquire_map(nimgs, rebin=rebinfact)



from splattsw import acceleration_analysis as acc
from splattsw.devices.webDAQ import WebDAQ

wbdq = WebDAQ()
wbdq.connect()

ffv = dm.mirrorModes
mode_id = 4
mode_amp = 5e-6
cmd = ffv[:,mode_id]*mode_amp

wbdq.start_schedule()
tn_buf = dm.sendBufferCommand()
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




from opticalib import analyzer as az
from opticalib.ground import zernike as zern
tn0='20251008_144454'
tn1='20251008_152943'
tn0='20251009_122248'
tn1='20251009_122904'

fl0 = osu.getFileList(tn0,key='20')
n0=len(fl0)
fl1 = osu.getFileList(tn1,key='20')
img00 = az.averageFrames(tn0,0,10)
img01 = az.averageFrames(tn0,n0-10,n0-1)
dd= zern.removeZernike(img00-img01,[1,2,3])
imshow(dd)
img10 = az.averageFrames(tn1,0,10)
#img01 = az.averageFrames(tn0,n0-10,n0-1)
dflat = zern.removeZernike(img10-img00,[1,2,3])


## thermo test
tn0='20251009_160205' #flat 50 frames, stable
tn1='20251009_160713' # rampup
tn2='20251009_162331' #flat
nframes = 50

img0 = zern.removeZernike(az.averageFrames(tn0,0,nframes),[1,2,3])
imshow(img0);title('Initial Flat: '+tn0);img0.std()
imgramp = zern.removeZernike(az.averageFrames(tn1,50,99)-az.averageFrames(tn1,0,50),[1,2,3])
imshow(imgramp);title('Diff after ramp '+tn1);
img1 = zern.removeZernike(az.averageFrames(tn2,0,nframes),[1,2,3])
imshow(img1);title('Flat after ramp '+tn2);img1.std()
dflat = img1-img0
imshow(dflat);title('Flat differences, '+tn2+'-'+tn0);dflat.std()





tn1 = '20251009_122904'# 727 frames ramp
img0 = zern.removeZernike(az.averageFrames(tn1,0,50),[1,2,3])
img1 = zern.removeZernike(az.averageFrames(tn1,650,700),[1,2,3])
imgramp=img1-img0
imshow(imgramp);title(tn1+' Diff after 5°C ramp,SfE='+str(int(imgramp.std()*1e9))+'nm');colorbar()
img0.std();img1.std()

tnim = '20251008_152616'
from opticalib.dmutils.flattening import Flattening
fl=Flattening(tnim)

img00 = az.modeRebinner(img0,4)
img11= az.modeRebinner(img1,4)
imf0=synthFlat(img00)
imf1=synthFlat(img11)
q0=img00+imf0
q1=img11+imf1
dd = q1-q0
w=20e-9
imshow(dd,vmin=-w,vmax=w);title('Flat diff. after ramp: SfE='+str(int(dd.std()*1e9))+'nm')

def synthFlat(img):
    fl.loadImage2Shape(img)
    fl.computeRecMat(3)
    rec=fl._recMat
    im = fl._loadIntCube()
    img = np.ma.masked_array(fl.shape2flat, mask=fl._getMasterMask())
    fc=-np.dot(img.compressed(), rec)
    print(fc)
    imf = img.copy()*0
    for i in range(19):
        imf += im[:,:,i]*fc[i]
    return imf

img1 = az.modeRebinner(img1,4)

imgstart = az.modeRebinner(img00, 4)
tnim = '20251008_152616'
from opticalib.dmutils.flattening import Flattening
fl=Flattening(tnim)
img2flat = az.modeRebinner(img10,4)
fl.loadImage2Shape(img2flat)
fl.computeRecMat(nmodes2discard)
deltacmd = fl.computeFlatCmd(nmodes2flat)

im=fl._loadIntCube()

img = np.ma.masked_array(fl.shape2flat, mask=fl._getMasterMask())
fc=-np.dot(img.compressed(), rec)
imf = img*0
for i in range(19):
    imf += im[:,:,i]*fc[i]


