#memo to run the SPLATT measurements in the laboratory
#importing
from SPLATT.splattsw import splatt4dmeas as comm4d
from SPLATT.splattsw import splatt_analysis as sp
from SPLATT.splattsw import splatt_acq as acq

from importlib  import reload #for reload
#power on and connect the RedPitaya
#start the redpitaya server
#systemctl start redpitaya_scpi 
from SPLATT.splattsw.devices.webDAQ import WebDAQ as wdq  
from SPLATT.splattsw.devices import redpitaya as rp
from M4.m4.mOTT_analysis import timehistory as th

webdaq= wdq()
webdaq.connect()
# end of initialization commands
#see cmdlog.py for the commands to be used

acq.acq(30,250)

freqM=np.array([20,240])
ampPI = 2
rp.splatt_trigger(freqM[0], freqM[1], ampPI)

webdaq.start_schedule()

#interferometer acquisition
tn = comm4d.capture(100)
#interferometer processing
comm4d.produce(tn)

wdfile= 'SPLATT_2021-11-24T17-12-59-351.wdd'

tn='20211020_165346'  #for testing
#interferometer analysis
tt = sp.tiltvec(tn)
tttimevec = sp.tttimevec(tn)
tt = tt[sp.tiltselect,:]
ttspe, ttfreq = sp.tt_spectrum(tt)
sp.ttplot(ttfreq, ttspe[sp.tiltselect,:])
sp.ttplot(ttfreq, ttspe)

tt = sp.tiltvec(tn)
tts, ttf = sp.tt_spectrum(tt)
tty = tts[sp.tiltselect,:]
ttpeak, ttpeakfreq = sp.find_peak(tty, freq=ttf, bound=bound)

f,v=sp.comb_analysis(f1, spe1[0,:],10)  #10 is the comb spacing

#accelerometer analysis
#read file
q = sp.openfile(wdfile)
#spectrum
spe, f = sp.acc_spectrum(v)
spe, f = sp.acc_spectrum(wdfile)
p.accplot(f, spe[0,:])


bound = np.array([18,25])
vobb = spe[0,:]
peak, peakfreq = sp.find_peak(vobb, freq=f, bound=bound)
wdfile=sp.lastwdfile()
sp.peak_analysis(tn, wdfile)


def peak_analysis(tn, wdfile, bound=np.array([18,25]):
    #bound = np.array([18,25])
    tt = sp.tiltvec(tn)
    tts, ttf = sp.tt_spectrum(tt)
    tty = tts[sp.tiltselect,:]
    ttpeak, ttpeakfreq = sp.find_peak(tty, freq=ttf, bound=bound)
    spe, f = sp.acc_spectrum(wdfile)
    peakOBB, peakfreqOBB = sp.find_peak(spe[0,:], freq=f, bound=bound)
    peakStand, peakfreqStand = sp.find_peak(spe[1,:], freq=f, bound=bound)
    peakOBBint = sp.acc_integrate(peakOBB, peakfreqOBB)
    peakStandint=sp.acc_integrate(peakStand, peakfreqStand)

    print('TT peak [m],   TT Freq [Hz]',)
    print(ttpeak, ttpeakfreq)

    print('OBB peak [m],  OBB Freq [Hz]')
    print(peakOBBint,peakfreqOBB)
    print('Stand peak [m],  Stand Freq [Hz]')
    print(peakStandint,peakfreqStand)

end
