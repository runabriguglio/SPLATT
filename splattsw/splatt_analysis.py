import os
import glob
import numpy as np
#import jdcal
from astropy.io import fits as pyfits
from SPLATT.splattsw.devices import webDAQ as wd
#from M4.m4.mOTT_analysis import timehistory as th
from M4.m4.mini_OTT import timehistory as th

from matplotlib.pyplot import *
rcParams['image.cmap'] = 'hot'
from SPLATT.splattsw import splatt_log as slog
from SPLATT.splattsw import splatt4dmeas as comm4d

#from m4.ground import read_data
#from m4.ground import zernike
#from m4.ground.read_data import InterferometerConverter
#ic = InterferometerConverter()
logfile = '/mnt/jumbo/SPLATT/'

basepathwebdaq= '/mnt/jumbo/SPLATT/WebDaqData/' #'/home/labot/ftp/'
freqwebdaq = 1651.6129 #Hz; minimum sampling frequency

ftpwebdacq = '/home/ftpuser/ftp/'

basepathbuffer = '/mnt/jumbo/SPLATT/Buffer/'
freqbuff = 118

basepath4d = '/mnt/jumbo/SPLATT/OPTData/'
freq4d = 122.64 #80.62
tn='tn0'
z2fit = np.array([1,2,3])
tiltselect = 1 #TBC

def wdsync():
    os.system('rsync -av '+ftpwebdacq+' '+basepathwebdaq)

def get_freq4d(tn):
    w = comm4d.read_acqdata(tn)
    f4d = w[0]    
    return f4d

def accplot(v, f):
    #plot(f, v)
    bar(f,v)
    xlabel('Frequency [Hz]')
    ylabel('Acceleration [g]')

def ttplot(v,f):
    bar(f,v)
    xlabel('Frequency [Hz]')
    ylabel('Tilt amplitude [m]')


def openfile(name):
    name = basepathwebdaq+name
    data = wd.openwdd(name)
    return data

def lastwdfile():
    #fl = th.fileList(basepathwebdaq)
    #fold1 = fold+'/'+tn+addfold   #to be re-checked at OTT!! 
    fl = sorted(glob.glob(basepathwebdaq+'SPLATT_Test*'))#, key=lambda x: int(os.path.basepathwebdaq(x).split('/')[-1].split('.')[0]))

    ffl = fl[-1]
    fname = ffl.split('/')[-1]
    print(fname)
    return fname


def acc_spectrum(v):
    mytype = type(v)
    if mytype is str:
        v1 = openfile(v)

    if mytype is np.array:
        v1=v
    
    if mytype is np.ndarray:
        v1=v

    spe, f = th.spectrum(v1, dt=1/freqwebdaq)
    return spe, f

def fileList(tn):
    flist=th.fileList(tn, fold=basepath4d)
    return(flist)

def tiltvec(tn):
    
    flist=th.fileList(tn, fold=basepath4d)
    tt = th.zernikePlot(flist, modes=z2fit)
    
    return tt

def tttimevec(tn):
    flist=th.fileList(tn, fold=basepath4d)
    tvec = np.arange(len(flist))/freq4d

    return tvec

def acc_timevec(v):
    q=openfile(v)
    npo=len(q[0,:])
    tve = np.arange(npo)/freqwebdaq
    return tve

def tt_spectrum(v, freq=None):
    mytype = type(v)
    if mytype is list:
        v1 = openfile(v)
        flist=th.fileList(tn, fold=basepath4d)
        v1 = th.zernikePlot(flist, modes=z2fit)

    if mytype is np.ndarray:
        v1=v

    mytype = type(freq)
    if mytype is str:
        freqsp = get_freq4d(freq)
        freq = freqsp
    if mytype is np.ndarray or mytype is np.int:
        freqsp=freq

    if freq is None: 
        freqsp=freq4d
    if freqsp == 0:
        freqsp = freq4d
    print(freqsp)
    spe, f = th.spectrum(v1, dt=1/freqsp)
    return spe, f

def wf_list(tn):
    flist=th.fileList(tn, fold=basepath4d)
    ww = []
    nn = np.size(flist)
    for i in np.arange(nn):
        img = th.frame(i,flist)
        ww.append(np.std(img))
    ww = np.asarray(ww)
    return ww

def wf_spectrum(tn, freq=None):

    wf = wf_list(tn) 

    mytype = type(freq)
    if mytype is str:
        freqsp = get_freq4d(freq)
        freq = freqsp
    if mytype is np.ndarray or mytype is np.int:
        freqsp=freq

    if freq is None:
        freqsp=freq4d
    if freqsp == 0:
        freqsp = freq4d
    print(freqsp)
    spe, f = th.spectrum(wf, dt=1/freqsp)
    return spe, f




def timevec():

    return tv

    
def runningMean(vec, npoints):
    
    return np.convolve(vec, np.ones(npoints), 'valid') / npoints       

def acc_integrate(acc, freq):
    #acc is the acceleration value in g
    #print('Is the acceleration provided as [g]?')
    acc = acc * 9.81 #converts to m/s2
    amp = acc/(4*np.pi**2*freq**2)
    return amp 


def find_peak(v, freq=None, bound=None):
    #requires a monodimensional array
    idf = range(len(v))
    if freq is not None and bound is not None:
        idf =     np.where(np.logical_and(freq>= bound[0], freq<=bound[1]))
    peak = max(v[idf])
    peakid = np.argmax(v[idf])
    peakfreq = 0
    if freq is not None:
        idf1 = idf[0]
        peakfreq = freq[idf1[peakid]]

    return peak, peakfreq

def peak_phase(v, freq, peakFreq):
    #requires a monodimensional array
    ph = phase_spectrum(v)
    idp = np.where(freq == peakFreq)
    pphase = ph[idp]
    pphase = pphase[0]
    return pphase

def phase_spectrum(v):
    sp = np.fft.fft(v)
    nn = int(np.size(sp)/2)
    sp = sp[0:nn+1]
    ph = np.angle(sp)
    return ph

def acc_peak_analysis(wdfile, bound=None, allaccel=None):
    wds         = openfile(wdfile)
    spe, f      = acc_spectrum(wds)
    #accphase    = phase_spectrum(wds)
    peakOBB, peakfreqOBB = find_peak(spe[0,:], freq=f, bound=bound)
    peakStand, peakfreqStand = find_peak(spe[1,:], freq=f, bound=bound)
    peakOBBint  = acc_integrate(peakOBB, peakfreqOBB)
    peakStandint=acc_integrate(peakStand, peakfreqStand)
    phaseOBB    = peak_phase(wds[0,:], f, peakfreqOBB)/np.pi*180 
    phaseStand  = peak_phase(wds[1,:], f, peakfreqStand)/np.pi*180
    phaseInput  = peak_phase(wds[2,:], f, peakfreqOBB)/np.pi*180
    if allaccel is not None:
        peakAcc2, peakfreqAcc2 = find_peak(spe[2,:], freq=f, bound=bound)
        peakAcc3, peakfreqAcc3 = find_peak(spe[3,:], freq=f, bound=bound)
        peakAcc2int  = acc_integrate(peakAcc2, peakfreqAcc2)
        peakAcc3int=acc_integrate(peakAcc3, peakfreqAcc3)
        peakAcc2int = peakAcc2int*1e6        
        peakAcc3int = peakAcc3int*1e6


    peakOBBint = peakOBBint*1e6
    peakStandint = peakStandint*1e6
    print('OBB peak [um],       OBB Freq [Hz]')
    print(peakOBBint,peakfreqOBB)
    print('Stand peak [um],     Stand Freq [Hz]')
    print(peakStandint,peakfreqStand)
    print('Input phase, OBB phase,   Stand phase')
    print(phaseInput, phaseOBB, phaseStand)
    
    plot(wds[0,:]/np.std(wds[0,:]))
    plot(wds[1,:]/np.std(wds[1,:]))
    plot(wds[2,:]/np.std(wds[2,:]))
    if allaccel is not None:
        return peakOBBint, peakStandint, peakAcc2int, peakAcc3int
    else:
        return peakOBBint, peakStandint, phaseInput, phaseOBB, phaseStand    

def opt_peak_analysis(tn, bound=None, sumtilt=None, wfe=None):
    if wfe is not None:
        tty = wf_list(tn)
    else:
        tt = tiltvec(tn)
    
    tts, ttf = tt_spectrum(tt, tn)
    if wfe is not None:
        tty = tts[themode,:]
    else:
        tty = tts
    if sumtilt is not None:
        tty = np.sqrt(tts[1,:]**2+tts[2,:]**2)
    ttpeak, ttpeakfreq = find_peak(tty, freq=ttf, bound=bound)
    phaseTT    = peak_phase(tt[tiltselect,:], ttf, ttpeakfreq)/np.pi*180
    ttpeak = ttpeak*1e6
    if sumtilt is not None:
        print('Summing X,Y tilt')
    print('TT peak [um],            TT Freq [Hz]')
    print(ttpeak, ttpeakfreq)
    return ttpeak, phaseTT

def peak_analysis(tn, wdfile=None, bound=None, sumtilt=None, wfe=None):
    tt_peak,tt_phase = opt_peak_analysis(tn, bound=bound, sumtilt=sumtilt, wfe=wfe)
    if wdfile is not None:
        accPeakOBB,accPeakStand, phIn, phOBB, phStand = acc_peak_analysis(wdfile, bound=bound)    
    else:
         accPeakOBB, accPeakStand,  phIn, phOBB, phStand = 0,0,0,0,1
    accInfo = comm4d.read_acqdata(tn)
    slog.log_data([tn, wdfile],[accInfo[0],accInfo[1],accInfo[2], accPeakOBB, accPeakStand, tt_peak, phIn, phOBB, phStand, tt_phase ])



#obsolete
def accpeak_analysis(wdfile, bound=np.array([4,300])):

    if wdfile is not None:
        spe, f = acc_spectrum(wdfile)
        peakOBB, peakfreqOBB = find_peak(spe[0,:], freq=f, bound=bound)
        peakStand, peakfreqStand = find_peak(spe[1,:], freq=f, bound=bound)
        peakOBBint = acc_integrate(peakOBB, peakfreqOBB)
        peakStandint=acc_integrate(peakStand, peakfreqStand)
        acc4 = 0
        if acc4 ==1: 
            peakSwing, peakfreqSwing = find_peak(spe[2,:], freq=f, bound=bound)
            peakBench, peakfreqBench = find_peak(spe[3,:], freq=f, bound=bound)
            peakSwingint = acc_integrate(peakSwing, peakfreqSwing)
            peakBenchint=acc_integrate(peakBench, peakfreqBench)

        print('OBB peak [m],       OBB Freq [Hz]')
        print(peakOBBint,peakfreqOBB)
        print('Stand peak [m],     Stand Freq [Hz]')
        print(peakStandint,peakfreqStand)
        if acc4 ==1:
            print('SwingArm peak [m],       SwingArm Freq [Hz]')
            print(peakSwingint,peakfreqSwing)
            print('Bench peak [m],     Bench Freq [Hz]')
            print(peakBenchint,peakfreqBench)
        #clf()
        #plot(f, spe[0,:])
        #plot(f, spe[1,:])
        #plot(f, spe[2,:])
        #plot(f, spe[3,:])
        return peakOBBint,peakStandint

def plotwaves(v):
    v1 = v[id1,:]
    v2 = v[id2,:]
    v = v/np.std(v1)
    
    plot(v1)
    plot(v2)



def readbuffer(tn):
    fname = basepathbuffer+tn
    hdu = pyfits.open(fname)
    data = hdu[0].data
    return data


def comp_freq2(f,v, threshold=None):
    dv= v[1:]-v[0:-1]
    
    if threshold is None:
        mythreshold = min(dv)/20
    else:
        mythreshold = threshold

    idf = where(dv < mythreshold)
    idf = idf[0]
    idf = idf-1
    return f[idf] 

def comb_freq(fsign, fref):
    comb_spacing = fsign
    fmax = np.max(fref)
    maxcomb=comb_spacing*(np.floor(fmax/comb_spacing))
    fout = np.arange(0+comb_spacing, fmax,comb_spacing)
    return fout

def extract_comb_id(f, comb ):
    theshold = f[1]-f[0]
    ncomb = np.size(comb)
    idout = []
    for i in comb:
        id1 = np.argmin(np.abs(f-i))
        idout.append(id1)

    
    return idout

def extract_comb(fid, v, freqbin=1):
    if freqbin is None:
        vout = v[fid]

    else:
        nfid = np.size(fid)
        vout = [] 
        fidout=[]
        for j in fid:
            #mv = []
            #for i in np.arange(-freqbin,freqbin+1):
            #    mv.append(v[j+i])
            #vout.append(np.mean(mv))
            v2eval = v[j+np.arange(-freqbin,freqbin+1)]
            vout.append(np.max(v2eval))
            fid2 = np.argmax(v2eval)
            fidout.append(j+fid2)
            #print(j,(v[j+np.arange(-freqbin,freqbin+1)]))
            
    vout = np.asarray(vout)
    fidout=np.asarray(fidout)
    return fidout, vout

def comb_analysis(f, v, fcomb,freqbin=1):
    
    comb= comb_freq(fcomb,f)   
    fid = extract_comb_id(f, comb )
    fout, vout =    extract_comb(fid, v, freqbin=5)
    return f[fout],vout
    
def sweep_analysis(tn, reb=None, tiltsum=None, wfe=None):
    if wfe is not None:
        v, f = wf_spectrum(tn, tn)
    else:
        tt = tiltvec(tn)
        spe, f = tt_spectrum(tt, tn)
        v = spe[1,:]
        if tiltsum is not None:
            t0=tiltsum[0]-1
            t1=tiltsum[1]-1
            v = np.sqrt(spe[t0,:]**2 + spe[t1,:]**2)

    if reb is not None:
        f = th.runningMean(f, reb)
        v = th.runningMean(v, reb)
    return f,v

def runningMean(v, npo):
    v1 = th.runningMean(v, npo)
    return v1
#def print_data(tn=tn, a=a, b=b):
#    print('%s\t%f\t%f' %(tn, a,b))
   

def frames_spectrum(tn):
    fl = th.fileList(tn)
    nframe= size(fl) 


def signal_unwrap(x, thr=632e-9/2, phase = 632e-9):
    v = x-x[0]
    npx = np.size(v)
    for i in np.arange(1,npx):
        dv = v[i]-v[i-1]
        if dv > thr:
            v[i] =v[i]-np.abs(phase * (round(dv/phase)))
        if dv < -thr:
            v[i] =v[i]+np.abs(phase * (round(dv/phase)))

    return v

def sweep_analysis_sequence(tnlist):
    lab = []
    ntn = np.size(tnlist)
    #cols = ['k','y','b','c','g','r','m','y']
    fl=th.fileList(tnlist[0], fold=basepath4d)
    npo=int(np.size(fl)/2+1)
    ff = np.zeros(npo)
    vv = np.zeros([npo,ntn])
    m2fit = np.array([2,3])
    for i in np.arange(ntn):
        f,v = sweep_analysis(tnlist[i],tiltsum=m2fit)
        #lab_i = ("Gap: %s" %gaplist[i])
        #lab.append(lab_i)
        ff=f
        vv[:,i] = v

    #fname = '/mnt/jumbo/SPLATT/OPTData/20220202_sweep_freq.fits'
    #pyfits.writeto(fname,ff)
    
    #for i in range(ntn-1):
    #    plot(ff,vv[:,i], label=lab[i])

    legend()
    return ff, vv

