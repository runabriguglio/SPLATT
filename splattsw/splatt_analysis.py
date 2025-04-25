import numpy as np
from aoptics import analyzer as th
from matplotlib.pyplot import *
z2fit = [1,2,3]
tiltselect = 1 #confirmed with PhaseCam4020 ex refurb 4030 and 4DFocus


def plot_ttspectra(spe,f, tn=None):
    figure()
    plot(f, spe[1,:])
    plot(f,spe[2,:])
    xlabel('Freq [Hz]')
    ylabel('Amp Spectrum [nm]')
    legend(['Tilt Y','Tilt X'])
    if tn is not None:
        title(tn)

def plot_multispec(spvec,f):
    ntn = len(tnlist)
    figure()
    for i in range(ntn):
        spe = spvec[i]
        plot(f, spe[1,:])
        #plot(f, spe[2,:])
    xlabel('Freq [Hz]')
    ylabel('Amp Spectrum [nm]')
    #legend(tnlist)



def tt_spectrum(tn, tt = None):
    flist=fileList(tn)
    freq4d = th.osu.getFrameRate(tn)
    if tt is None:
        print('Computing Z vec')
        tt = tiltvec(tn)
#aggiungere unwrap
    spe, f = th.spectrum(tt, dt=1/freq4d)
    return spe ,f



def tiltvec(tn, unwrap=True):
    flist=fileList(tn)
    nf = len(flist)
    nm = len(z2fit)
    tt = th.zernikePlot(flist, modes=z2fit)
    if unwrap == True:
        pt = tt[0,:]
        p0 = signal_unwrap(pt)
        tt1 = np.zeros([nm+1,nf])
        tt1[0:nm,:] = tt
        tt1[nm,:] = p0
    else:
        tt1 = tt
    return tt1

def timevec(tn):
    flist=fileList(tn)
    freq4d = th.osu.getFrameRate(tn)
    tvec = np.arange(len(flist))/freq4d
    return tvec

def fileList(tn):
    flist=th.osu.getFileList(tn, key=".4D")
    return flist

def averageFrames(tn,first=None,last=None):
    flist=th.osu.getFileList(tn, key=".4D")
    if first is None and last is None:
        cc = th.createCube(flist)
    else:
        cc = th.createCube(flist[first:last])
    aveimg = np.ma.average(cc,2)
    return aveimg


def runningMean(vec, npoints):
    v1 = th.runningMean(vec, npoints)
    return v1

def signal_unwrap(x, thr=632e-9/4, phase = 632e-9/2):
    v = x-x[0]
    npx = np.size(v)
    for i in np.arange(1,npx):
        dv = v[i]-v[i-1]
        if dv > thr:
            v[i] =v[i]-np.abs(phase * (round(dv/phase)))
        if dv < -thr:
            v[i] =v[i]+np.abs(phase * (round(dv/phase)))
    return v


