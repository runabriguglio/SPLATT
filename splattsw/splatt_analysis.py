import numpy as np
from aoptics import analyzer as th
from matplotlib.pyplot import *
import json, struct
z2fit = [1,2,3]
tiltselect = 1 #confirmed with PhaseCam4020 ex refurb 4030 and 4DFocus
basepathwebdaq= '/mnt/jumbo/SPLATT/WebDaqData/' #'/home/labot/ftp/'
freqwebdaq = 1651.6129 #Hz; minimum sampling f



def plot_ttspectra(spe,f, tn=None):
    figure()
    plot(f, spe[1,:])
    plot(f,spe[2,:])
    xlabel('Freq [Hz]')
    ylabel('Amp Spectrum [nm]')
    legend(['Tilt Y','Tilt X'])
    if tn is not None:
        title(tn)



def tt_spectrum(tn, tt = None):
    flist=fileList(tn)
    freq4d = th.osu.getFrameRate(tn)
    if tt is None:
        print('Computing Z vec')
        tt = tiltvec(tn)
#aggiungere unwrap
    spe, f = th.spectrum(tt, dt=1/freq4d)
    return spe ,f

def acc_spectrum(wdname):
    q = openwdfile(wdname)
    spe, f = th.spectrum(q, dt=1/freqwebdaq)
    return spe, f

def mech_spectrum(wdname):
    spe, f = acc_spectrum(wdname)
    npacc = np.shape(spe)
    spei = np.zeros(npacc)
    for i in range(npacc[0]):
        v = spe[i,:]
        s = acc_integrate(v,f)
        spei[i,:] =s
    spei = np.array(spei)
    return spei, f

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



def runningMean(v, npo):
    v1 = th.runningMean(v, npo)
    #v1 = np.convolve(v,np.ones(npo)/npo)
    return v1

def signal_unwrap(x, thr=632e-9/4, phase = 632e-9/2):
    # v = np.unwrap(x-x[0], period = 632e-9/2)
    v = x-x[0]
    npx = np.size(v)
    for i in np.arange(1,npx):
        dv = v[i]-v[i-1]
        if dv > thr:
            v[i] =v[i]-np.abs(phase * (round(dv/phase)))
        if dv < -thr:
            v[i] =v[i]+np.abs(phase * (round(dv/phase)))
    return v


def openwdfile(thename):
    fname = basepathwebdaq+thename
    hdr = {}
    with open(fname,"rb") as wdf:
        rawData = wdf.read(564)
        hdr['version'] = int.from_bytes(rawData[0:4], 'little')
        hdr['size'] = int.from_bytes(rawData[4:8], 'little')
        hdr['nchannels']= int.from_bytes(rawData[8:12], 'little')
        hdr['scan_rate']= int.from_bytes(rawData[12:20], 'little')
        hdr['start_time']= int.from_bytes(rawData[20:28], 'little')
        hdr['timezone']= rawData[28:44].decode("utf-8")
        json_hdr_size =  int.from_bytes(rawData[560:564], 'little')
        jsonRaw = wdf.read(json_hdr_size)
        hdr['json_hdr']=json.loads(jsonRaw)
        ndata = hdr['json_hdr']['jobDescriptor']['acquisition']['stopTrigger']['sampleCount']
        data_it = struct.iter_unpack('<d', wdf.read(ndata*hdr['nchannels']*8)) #4 because double precision 64 bit\n",
        tmp = np.asarray(list(data_it), dtype='double')
        data=tmp.reshape(int(tmp.size/4), 4)
        data = data.T
    return data

def acc_integrate(acc, freq):
    #acc is the acceleration value in g
    #print('Is the acceleration provided as [g]?')
    acc = acc * 9.81 #converts to m/s2
    amp = acc/(4*np.pi**2*freq**2)
    return amp
