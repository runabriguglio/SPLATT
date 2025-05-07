import os, numpy as np
from opticalib import analyzer as th
from matplotlib.pyplot import *
import json, struct
from astropy.io import fits as pyfits
from matplotlib.pyplot import *



z2fit            = [1,2,3]
tiltselect       = 1 #confirmed with PhaseCam4020 ex refurb 4030 and 4DFocus
basepathwebdaq   = '/mnt/jumbo/SPLATT/WebDaqData/' #'/home/labot/ftp/'
freqwebdaq       = 1651.6129 #Hz; minimum sampling f
foldanalysisconf = '/home/labot/git/SPLATT/'
resultfold       = '/mnt/libero/SPLATTData/Results/'

#configureRC()

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
    spe, f = th.spectrum(tt, dt=1/freq4d)
    return spe ,f

def acc_spectrum(wdname):
    q = openwdfile(wdname)
    spe, f = th.spectrum(q, dt=1/freqwebdaq)
    return spe, f

def mech_spectrum(accvec):
    spe, f = th.spectrum(accvec, dt=1/freqwebdaq)
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

def zvec(tn, overwrite = False):
    resultpath = resultfold+tn+'/'
    if overwrite == True or os.path.exists(resultpath) == False:
        os.mkdir(resultpath)
        flist=fileList(tn)
        nfiles = nel(flist)
        imgcube = th.createCube(flist)
        zzvec = np.zeros(nzern, nfiles)
        rrvec = np.zeros(nfiles)
        for i in range(nfiles):
            cc, m = th.zern.zernikeFit(imgcube[i],nzern)
            zzvec[:,i] = cc
            zr = np.sum(cc**2)
            rrvec[i] = np.sqrt(imgcube[i].std()**2 - zr)
        pyfits.writeto(resultpath+'-ZernVec.fits',zzvec)
        pyfits.writeto(resultpath+'-stdVec.fits', rrvec)
    else:
        h0 = pyfits.open(resultpath+'-ZernVec.fits')
        zzvec = h0[0].data
        h0.close()
        h0 = pyfits.open(resultpath+'-stdVec.fits')
        rrvec = h0[0].data
        h0.close()
    return zzvec, rrvec




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


def openaccfile(tn):
    accfile = basepathwebdaq+ tn+'.fits'
    h0 = pyfits.open(accfile)
    accelerations = h0[0].data
    h0.close()
    return accelerations


def acc_integrate(acc, freq):
    #acc is the acceleration value in g
    #print('Is the acceleration provided as [g]?')
    acc = acc * 9.81 #converts to m/s2
    amp = acc/(4*np.pi**2*freq**2)
    return amp


def read_analysisconf(tn):
    import configparser
    config=configparser.ConfigParser()
    config.read(foldanalysisconf+tn+'.meas')
    dataset = config['DATASET']
    meas    = config['MEASUREMENT']
    analysis= config['ANALYSIS']
    
    
    return dataset, meas, analysis

def sweep_analysis(tn, nrunnmean: int = 5):
    '''
    speczrip: 
    specsrip
    spec
    '''
    savefold = resultfold+tn+'/'
    if os.path.exists(savefold) is False:
        os.mkdir(savefold)
    tninfo,measinfo, analysisinfo = read_analysisconf(tn)
    '''
    content of theconf:
        tnrip
        tnset: list of tracknums. each tracknum includes optical data, accelerometer data, buffer and test config
        tnlabel: list (same len as tnset) of identifier to be associated with tnset tracknums. is to be used as labels for the plot
    '''
    tnrip   = tninfo['tnrip']
    tnset   = tninfo['tnset']
    tnlabel = tninfo['tnlabel']
    
    ### acceleration analysis
    accv       = openaccfile(tnrip)
    acctimevec = np.arange(len(accv[0]))/freqwebdaq

    mspec, mfreq = mech_spectrum(accv)

    #utils
    accxlims = json.loads(analysisinfo['acclim'])
    acclegend = json.loads(tninfo['acclabel'])
    figure(figsize=(12,6));          suptitle(tn+' Accelerometer data')
    subplot(1,2,1);    title('Acceler. time series')
    for i in accv:
        plot(acctimevec, i,'.')
    xlabel('Time [s]');    ylabel('Acceleration [g]');    legend(acclegend)

    subplot(1,2,2);     title('Mechanical Oscillation spectrum')
    for i in mspec:
        plot(mfreq, i,'.')
    xlabel('Freq [Hz]'); ylabel(' Amplitude [m]'); legend(acclegend); yscale('log');xlim(accxlims[0],accxlims[1])
    savefig(savefold+tn+'-AccPlots.png')
    ###    end of acc data    ############
    pass
    ## analysis of optical data ###
    ntn        = len(tnset)
    freq4d     = th.osu.getFrameRate(tn)
    zv_rip, st_rip = zvec(tnrip)
    speczrip, f = th.spectrum(zv_rip, dt=1/freq4d) 
    specsrip, f = th.spectrum(st_rip, dt=1/freq4d)
    frunn       = th.runningmean(f, nrunnmean)
    speczrip    = th.runningmean(speczrip, nrunnmean)
    specsrip    = th.runningmean(specsrip, nrunnmean)
    zv       = []
    sv       = []
    speczset = []
    specsset = []
    spenz    = []
    spens    = []
    for i in tnset:
        zvi, svi = zvec(i)
        zv.append(zvi)
        sv.append(svi)
        specz, f = th.spectrum(zvi, dt=1/freq4d)
        specz    = th.runningmean(specz, nrunnmean)
        specs, f = th.spectrum(svi, dt=1/freq4d)
        specs    = th.runningmean(specs, nrunnmean)
        speczset.append(specz)
        specsset.append(specs)
        spenz.append(specz/speczrip)
        spens.append(specs/speczrip)

    ## plot of single spectra
    ### Excited Tilt
    figure(); suptitle(tn+' Optical SurfError')
    subplots_adjust(bottom=0.1, left=0.1,right=0.95, top=0.9,wspace=0.3, hspace=0)
    subplot(121)
    plot(frunn, speczrip[tiltselect,:],'.')
    for i in range(tnset):
        thespec = speczset[i]
        plot(frunn, thespec[tiltselect,:],'.')
    title(tn+' Y Opt. tilt amp spectrum')
    xlabel('Freq [Hz]');    ylabel('Tilt amp. [m]');    legend(['RIP',tnlabel])
    yscale('log')

    ### Global Tilt
    #??? è corretto sommare in quadr. le ampiezze degli spettri???
    tiltrip = np.sqrt(speczrip[1,:]**2+ speczrip[2,:]**2)
    subplot(222)
    plot(frunn, tiltrip,'.')
    for i in range(tnset):
        thespec = speczset[i]
        thetilt = np.sqrt(thespec[1,:]**2+ thespec[2,:]**2)
        plot(frunn, thespec,'.')
    title(tn+' Opt. tilt amp spectrum')
    xlabel('Freq [Hz]');    ylabel('Tilt amp. [m]');    legend(['RIP',tnlabel])
    yscale('log')


    ### Astigmatism
    #??? è corretto sommare in quadr. le ampiezze degli spettri???
    asttrip = np.sqrt(speczrip[4,:]**2+ speczrip[5,:]**2)
    subplot(223)
    plot(frunn, tiltrip,'.')
    for i in range(tnset):
        thespec = speczset[i]
        theast = np.sqrt(thespec[4,:]**2+ thespec[5,:]**2)
        plot(frunn, theast,'.')
    title(tn+' Opt. Astigm amp spectrum')
    xlabel('Freq [Hz]');    ylabel('Astigm amp. [m]');    legend(['RIP',tnlabel])
    yscale('log')


    ### Global RMS
    #??? è corretto sommare in quadr. le ampiezze degli spettri???
    subplot(224)
    plot(frunn, specsrip,'.')
    for i in range(tnset):
        thespec = specsset[i]
        plot(frunn, thespec,'.')
    title(tn+' Res. SFE spectrum')
    xlabel('Freq [Hz]');    ylabel('SFE [m]');    legend(['RIP',tnlabel])
    yscale('log')

    ##Normalized dataset
    ### Excited Tilt
    figure();  suptitle(tn+' Optical spectra, norm.')
    subplot(121); title('Y Opt. tilt amp spectrum')
    plot(frunn, np.ones(len(frunn)),'.')
    for i in range(tnset):
        thespec = specnz[i]
        plot(frunn, thespec[tiltselect,:],'.')
    xlabel('Freq [Hz]');    ylabel('Tilt norm. amp. [m]');    legend(['REF',tnlabel])
    yscale('log')


    ### Normalized RMS
    #??? è corretto sommare in quadr. le ampiezze degli spettri???
    subplot(122);  title('Res. SFE spectrum')
    plot(frunn, np.ones(len(frunn)),'.')
    for i in range(tnset):
        thespec = specns[i]
        plot(frunn, thespec,'.')
    xlabel('Freq [Hz]');    ylabel('SFE-norm. [m]');    legend(['REF',tnlabel])
    yscale('log')


def singlefreq_analysis(tn):
    pass


def configureRC():
    ks = ['top','bottom','left','right','wspace','hspace']
    val = [0.9,0.1,0.1,0.95,0.2,0.1]
    npar = len(ks)
    for i in range(npar):
        nname = 'figure.subplot.'+ks[i]
        rcParams[nname] = val[i]

configureRC()
