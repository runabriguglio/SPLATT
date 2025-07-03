import os, numpy as np
from opticalib import analyzer as th
from matplotlib.pyplot import *
import json, struct
from astropy.io import fits as pyfits
from matplotlib.pyplot import *

#this part shall be moved to acc_analysis
pmat = np.array([[1,1],[1,-1]])
p1mat = np.linalg.inv(pmat)

z2fit            = [1,2,3]
tiltselect       = 1 #confirmed with PhaseCam4020 ex refurb 4030 and 4DFocus
basepathwebdaq   = '/mnt/jumbo/SPLATT/WebDaqData/' #'/home/labot/ftp/'
freqwebdaq       = 1651.6129 #Hz; minimum sampling f
foldanalysisconf = '/home/labot/git/SPLATT/'
resultfold       = '/mnt/libero/SPLATTData/Results/'
testconfigpath   = '/mnt/libero/SPLATTData/TestConfig/'
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
    freq4d = read_4dfreq(tn) #th.osu.getFrameRate(tn)
    if tt is None:
        print('Computing Z vec')
        tt, sv = zvec(tn)
    spe, f = th.spectrum(tt, dt=1/freq4d)
    return spe ,f

def acc_spectrum(tn):
    thetype = type(tn)
    if thetype == str:
        q = openaccfile(tn)
    else:
        q = tn.copy()
    spe, f = th.spectrum(q, dt=1/freqwebdaq)
    return spe, f

def acc2pt(q):
    qq = q.copy()
    qp = p1mat @ q[0:2,:]
    qq[0:2,:] = qp
    return qq

def mech_spectrum(tn,transform = True):
    q = openaccfile(tn)
    if transform is True:
       q = acc2pt(q) 
    spe, f = acc_spectrum(q)
    npacc = np.shape(spe)
    spei = np.zeros(npacc)
    for i in range(npacc[0]):
        v = spe[i,:]
        s = acc_integrate(v,f)
        spei[i,:] =s
    spei = np.array(spei)
    return spei, f

def plot_mechspectrum(tn, f0=None):
    figure()
    acclabelT=["Pist","TiltY","Stand","ElevArm"]
    spm, fm = mech_spectrum(tn,transform = True)
    for i in spm:
        plot(fm, i,'.')
    yscale('log')
    xlim(0,120);grid('on');xlabel('Freq [Hz]');ylabel('Displac. amplitude spectrum [m Peak]');title('Accelerom. data: '+tn)
    legend(acclabelT)
    dp=[]
    ds = []
    mybound = [2,120]
    if f0 is not None:
        mybound =[f0-1,f0+1]
        for i in spm:
            pm, fmm = find_peak(i,freq=fm, bound=mybound,integrate=False)
            print('Found peak at Freq:'+str(fmm)+' Hz')
            dp.append(pm)
        dp=np.array(dp)
        for i in range(4):
            ds.append(acclabelT[i]+': '+f'{dp[i]*1e9:.0f}'+'nm')

        legend(ds)
    return dp




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
    print('Processing dataset: '+tn)
    print('Double check the StD computation')
    resultpath = resultfold+tn
    z2fit = [1,2,3,4,5,6]
    nzern = len(z2fit)
    if (overwrite == True) or os.path.exists(resultpath+'-ZernVec.fits') == False:
        flist=fileList(tn)
        nfiles = len(flist)
        imgcube = th.createCube(flist)
        zzvec = np.zeros([nzern+1, nfiles])
        rrvec = np.zeros(nfiles)
        for i in range(nfiles):
            cc, m = th.zern.zernikeFit(imgcube[:,:,i],z2fit)
            zzvec[0:nzern,i] = cc
            zr = np.sum(cc[1:]**2)
            rrvec[i] = np.sqrt(imgcube[:,:,i].std()**2 - zr)
        zzvec[nzern,:]= signal_unwrap(zzvec[0,:])
        pyfits.writeto(resultpath+'-ZernVec.fits',zzvec, overwrite=True)
        pyfits.writeto(resultpath+'-stdVec.fits', rrvec, overwrite=True)
    else:
        h0 = pyfits.open(resultpath+'-ZernVec.fits')
        zzvec = h0[0].data
        h0.close()
        h0 = pyfits.open(resultpath+'-stdVec.fits')
        rrvec = h0[0].data
        h0.close()
    return zzvec, rrvec


def find_peak(v, freq=None, bound=None, integrate=False):
    #requires a monodimensional array
    idf = range(len(v))
    if freq is not None and bound is not None:
        idf =     np.where(np.logical_and(freq>= bound[0], freq<=bound[1]))
        #idf = idf[0]
        #print(len(idf[0]))
    peak = max(v[idf])
    #peakfreq = np.argmax(v[idf])
    #print(peakfreq)
    #print(idf[peakfreq])
    #print('Found peak at Freq: ')
    #print(freq[idf[peakfreq]])
    if integrate == True:
        peak = np.sum(v[idf])
    peakid = np.argmax(v[idf])
    peakfreq = 0
    if freq is not None:
        idf1 = idf[0]
        peakfreq = freq[idf1[peakid]]

    return peak, peakfreq



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



def runningMean(v, npo,dim=2):
    nn = np.shape(v)
    if dim == 2:
        vout = []
        for i in range(nn[dim-2]):
            vout.append(th.runningMean(v[i,:],npo))
    if dim == 1:
        vout = []
        for i in range(nn[dim-1]):
            vout.append(th.runningMean(v[:,i],npo))
    
    if dim == 0:
        vout = th.runningMean(v,npo)
    vout = np.array(vout)
    #v1 = np.convolve(v,np.ones(npo)/npo)
    return vout

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
    print(tn)
    accfile = basepathwebdaq+ tn+'.fits'
    print('Opening file:'+accfile)
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

def read_4dfreq(tn):
    fname = os.path.join(testconfigpath,tn,'frequency4D.fits')
    if os.path.exists(fname) == 1:
        print('Reading 4D Freq from log file')
        freq = (pyfits.open(fname))[0].data
        freq = freq[0]
    else:
        freq = th.osu.getFrameRate(tn)
    print('4D Freq:'+str(freq))
    return freq

def read_analysisconf(tn):
    import configparser
    config=configparser.ConfigParser()
    config.read(foldanalysisconf+tn+'.meas')
    dataset = config['DATASET']
    meas    = config['MEASUREMENT']
    analysis= config['ANALYSIS']
    
    
    return dataset, meas, analysis

def compare_tn(tn0, tn1, meastype = 'sweep', freq=None, nbins=9):
    frunn0, specz0, specs0 = dataprocess(tn0,meastype, freq, nbins)
    frunn1, specz1, specs1 = dataprocess(tn1,meastype, freq, nbins)
    plot(frunn0,specz0[1,:],'.')
    plot(frunn1,specz1[1,:],'.')
    yscale('log'); xlim(20,120);grid()
    legend([tn0, tn1])



def dataprocess(tn, meastype = 'sweep', freq=None, nbins=1):
    if meastype == 'single' and freq == None:
        print('Single freq measurement, frequency is requested')
        return
    zv, sv = zvec(tn)
    freq4d = read_4dfreq(tn)#implementare la lettura del file da TestConfig
    spe, f = th.spectrum(zv,dt=1/freq4d)
    spes,f = th.spectrum(sv, dt=1/freq4d)
    if meastype == 'sweep':
        frunn    = runningMean(f, nbins, dim=0)
        specz    = runningMean(spe, nbins)
        specs    = runningMean(spes, nbins,dim=0)
        return frunn, specz, specs

    if meastype == 'single':
        peaks, pf = find_peak(spes,f,bound = [freq-nbins,freq+nbins])
        pvec = []
        for m in range(len(zv)):
            peakz, pf = find_peak(spe[m,:],f,bound = [freq-nbins,freq+nbins])
            pvec.append(peakz)
            print(pf)
        return freq, pvec, peaks


def global_analysis(tn):
    '''
    speczrip: 
    specsrip
    spec
    '''
    #savefold = resultfold+tn+'/'
    #if os.path.exists(savefold) is False:
    #    os.mkdir(savefold)
    tninfo,measinfo, analysisinfo = read_analysisconf(tn)
    '''
    content of theconf:
        tnrip
        tnset: list of tracknums. each tracknum includes optical data, accelerometer data, buffer and test config
        tnlabel: list (same len as tnset) of identifier to be associated with tnset tracknums. is to be used as labels for the plot
    '''
    meastype = measinfo['type']
    if meastype == 'single':
        ptsym = '-o'
        ptrsym= '-x'
    if meastype == 'sweep':
        ptsym = '.'
        ptrsym= 'x'

    tnrip   = json.loads(tninfo['tnrip'])
    #if meastype == "single":
    #    tnrip = tnrip[0]
    tnset   = json.loads(tninfo['tnset'])
    tnlabel = json.loads(tninfo['tnlabel'])
    optlabel = []
    optlabel.append( 'RIP')
    optlabel.extend(tnlabel)
    optxlims = json.loads(analysisinfo['optlim'])
    figtitle = tninfo['testname']
    
    ### acceleration analysis
    if measinfo['accdata'] == 'yes':
        print('Processing accel. data')
        dotransform = False
        if analysisinfo['proj_acc'] == 'yes':
            print('Projecting acc signals')
            dotransform = True
        print(tnrip)
        accv       = openaccfile(tnrip)
        acctimevec = np.arange(len(accv[0]))/freqwebdaq

        mspec, mfreq = mech_spectrum(tnrip, transform = dotransform)

    #utils
        accxlims = json.loads(analysisinfo['acclim'])
        acclegend = json.loads(tninfo['acclabel'])
        
        figure(10,figsize=(12,6));          suptitle(tn+' Accelerometer data')
        subplot(1,2,1);    title('Acceler. time series')
        for i in accv:
            plot(acctimevec, i,'.')
        xlabel('Time [s]');    ylabel('Acceleration [g]');    legend(acclegend)

        subplot(1,2,2);     title('Mechanical Oscillation spectrum')
        if analysisinfo['proj_acc'] == 'yes':
            acclegend = json.loads(tninfo['acclabel_proj'])
        for i in mspec:
            plot(mfreq, i,'.')
        xlabel('Freq [Hz]'); ylabel(' Amplitude [m]'); legend(acclegend); yscale('log');xlim(accxlims[0],accxlims[1]);grid()
        savefig(resultfold+tn+'-AccPlots.png')
    ###    end of acc data    ############
    
    ## analysis of optical data ###
    ntn        = len(tnset)

    zv       = []
    sv       = []
    speczset = []
    specsset = []
    spenz    = []
    spens    = []

    #analysis of RIP data
    if meastype == 'single':
        frunn = json.loads(measinfo['freqvector'])
        peakbins = json.loads(analysisinfo['peakpoints'])

        speczrip = []
        specsrip = []
        for i in range(len(tnrip)):
            fpeak,pvec, peaks = dataprocess(tnrip[i],meastype=meastype,freq=frunn[i],nbins=peakbins)
            speczrip.append(pvec)
            specsrip.append(peaks)
        speczrip = np.array(speczrip).T
        specsrip = np.array(specsrip).T

        for i in tnset:
            specz = []
            specs = []
            for m in range(len(i)):
                fpeak,pvec, peaks = dataprocess(i[m],meastype=meastype,freq=frunn[m],nbins=peakbins)
                specz.append(pvec)
                specs.append(peaks)
            specz = np.array(specz).T
            specs = np.array(specs).T

            speczset.append(specz)
            specsset.append(specs)
            spenz.append(specz/speczrip)
            spens.append(specs/specsrip)


    if meastype == 'sweep':
        nbins = json.loads(analysisinfo['runnmeanpoints'])
        print(nbins)
        frunn, speczrip, specsrip = dataprocess(tnrip,meastype=meastype,freq=None,nbins=nbins)

        for i in tnset:
            frunn, specz, specs = dataprocess(i,meastype=meastype,freq=None,nbins=nbins)
            speczset.append(specz)
            specsset.append(specs)
            spenz.append(specz/speczrip)
            spens.append(specs/specsrip)


        #original code
        '''
        zv_rip, st_rip = zvec(tnrip)
        freq4d     = th.osu.getFrameRate(tnrip)
        speczrip, f = th.spectrum(zv_rip, dt=1/freq4d) 
        specsrip, f = th.spectrum(st_rip, dt=1/freq4d)
        frunn       = runningMean(f, nrunnmean, dim=0)
        speczrip    = runningMean(speczrip, nrunnmean)
        specsrip    = runningMean(specsrip, nrunnmean,dim=0)

        zvi, svi = zvec(i)
        zv.append(zvi)
        sv.append(svi)
        specz, f = th.spectrum(zvi, dt=1/freq4d)
        specz    = runningMean(specz, nrunnmean)
        specs, f = th.spectrum(svi, dt=1/freq4d)
        specs    = runningMean(specs, nrunnmean,dim=0)
        speczset.append(specz)
        specsset.append(specs)
        spenz.append(specz/speczrip)
        spens.append(specs/specsrip)
        '''



    #resulting variables after this stage:
    #frunn: vector of frequencies
    # RIP: 
    #    speczrip (Zern spectrum)
    #    specsrip (stddev spectrum)
    # SET:
    #    speczset (Zern spectrum)
    #    specsset (stddev spectrum)
    
    ## plot of single spectra
    ### Excited Tilt
    figure(figsize=(10,5)); suptitle(figtitle+'\nOptical SurfaceError')
    subplots_adjust(bottom=0.15, left=0.1,right=0.95, top=0.85,wspace=0.3, hspace=0)
    subplot(121)
    plot(frunn, speczrip[tiltselect,:],ptsym)
    for i in range(ntn):
        thespec = speczset[i]
        plot(frunn, thespec[tiltselect,:],ptsym)
    title('Y Opt. tilt amp spectrum')
    xlabel('Freq [Hz]\nDataset: '+tn);    ylabel('Tilt amp. [m]');    legend(optlabel); xlim(optxlims[0],optxlims[1])
    yscale('log');grid()


    ### Global Tilt
    #??? è corretto sommare in quadr. le ampiezze degli spettri???
    #no. anche perchè la somma in quadr è definitia positiva e questo falserebbe il contenuto in frequenze
    #tiltrip = np.sqrt(speczrip[1,:]**2+ speczrip[2,:]**2)
    subplot(122)
    plot(frunn, speczrip[2,:],ptsym)
    for i in range(ntn):
        thespec = speczset[i]
        #thetilt = np.sqrt(thespec[1,:]**2+ thespec[2,:]**2)
        plot(frunn, thespec[2,:],ptsym)
    title('X Opt. tilt amp spectrum')
    xlabel('Freq [Hz]\nDataset: '+tn);    ylabel('Tilt amp. [m]');    legend(optlabel); xlim(optxlims[0],optxlims[1])
    yscale('log');grid()
    savefig(resultfold+tn+'-OptTilt.png')


    ## PLOT2
    ### Global RMS
    #??? è corretto sommare in quadr. le ampiezze degli spettri???
    figure(figsize=(15,5));suptitle(figtitle+'\nOptical Surface Error')
    subplots_adjust(bottom=0.15, left=0.1,right=0.95, top=0.85,wspace=0.3, hspace=0)
    subplot(131)
    plot(frunn, specsrip,ptsym)
    for i in range(ntn):
        thespec = specsset[i]
        plot(frunn, thespec,ptsym)
    title('Res. SFE spectrum')
    xlabel('Freq [Hz]\nDataset: '+tn);    ylabel('SFE [m]');    legend(optlabel); xlim(optxlims[0],optxlims[1])
    yscale('log');grid()

    ### Astigmatism
    #??? è corretto sommare in quadr. le ampiezze degli spettri???
    #asttrip = np.sqrt(speczrip[4,:]**2+ speczrip[5,:]**2)
    subplot(132)
    plot(frunn, speczrip[4,:],ptsym)
    for i in range(ntn):
        thespec = speczset[i]
        #theast = np.sqrt(thespec[4,:]**2+ thespec[5,:]**2)
        plot(frunn, thespec[4,:],ptsym)
    title('X Opt. Astigm amp spectrum')
    xlabel('Freq [Hz]\nDataset: '+tn);    ylabel('Astigm amp. [m]');    legend(optlabel); xlim(optxlims[0],optxlims[1])
    yscale('log');grid()

    subplot(133)
    plot(frunn, speczrip[5,:],ptsym)
    for i in range(ntn):
        thespec = speczset[i]
        #theast = np.sqrt(thespec[4,:]**2+ thespec[5,:]**2)
        plot(frunn, thespec[5,:],ptsym)
    title('Y Opt. Astigm amp spectrum')
    xlabel('Freq [Hz]\nDataset: '+tn);    ylabel('Astigm amp. [m]');    legend(optlabel); xlim(optxlims[0],optxlims[1])
    yscale('log');grid()
    savefig(resultfold+tn+'-OptAst-SfE.png')

    ##Normalized dataset
    optlabel[0] = 'REF'
    ### Excited Tilt
    figure(figsize=(9,5));  suptitle(figtitle+'\nOptical spectra, normalized to RIP');    subplots_adjust(bottom=0.15, left=0.1,right=0.95, top=0.85,wspace=0.3, hspace=0)

    subplot(121); title('Y Opt. tilt amp spectrum')
    plot(frunn, np.ones(len(frunn)),ptrsym)
    for i in range(ntn):
        thespec = spenz[i]
        plot(frunn, thespec[tiltselect,:],ptsym)
    xlabel('Freq [Hz]\nDataset: '+tn);    ylabel('Tilt norm. amp. [m]');    legend(optlabel); xlim(optxlims[0],optxlims[1])
    yscale('log');grid()


    ### Normalized RMS
    #??? è corretto sommare in quadr. le ampiezze degli spettri???
    subplot(122);  title('Res. SFE spectrum')
    plot(frunn, np.ones(len(frunn)),ptrsym)
    for i in range(ntn):
        thespec = spens[i]
        plot(frunn, thespec,ptsym)
    xlabel('Freq [Hz]\nDataset: '+tn);    ylabel('SFE-norm. [m]');    legend(optlabel); xlim(optxlims[0],optxlims[1])
    yscale('log');grid()
    savefig(resultfold+tn+'-OptNorm.png')

    return tninfo,measinfo, analysisinfo




def configureRC():
    ks = ['top','bottom','left','right','wspace','hspace']
    val = [0.9,0.1,0.1,0.95,0.2,0.1]
    npar = len(ks)
    for i in range(npar):
        nname = 'figure.subplot.'+ks[i]
        rcParams[nname] = val[i]

configureRC()
