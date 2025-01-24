
from splattsw import splatt4dmeas as comm4d
from splattsw import splatt_analysis as sp
from importlib  import reload #for reload
#power on and connect the RedPitaya
#start the redpitaya server
#systemctl start redpitaya_scpi 
from splattsw.devices.webDAQ import WebDAQ as wdq
from splattsw.devices import redpitaya as rp
from splattsw import splatt_log as slog
from matplotlib.pyplot import *

import os
import numpy as np
import time
webdaq= wdq()
webdaq.connect()
nframes = 2000
pauseStart = 10

def acq(freqPI, freq4d, ampPI=2, bound=None,process=None):
    nfr = nframes
    #rp.splatt_trigger(freqPI, freq4d, ampPI)
    rp.splatt_trigger2(freqPI, freq4d, ampPI)  # with pulse
    print('Command to RedPitaya: done!')
    print('Now waiting for transitory oscillations to damp')
    time.sleep(pauseStart)
    webdaq.start_schedule()
    
    tn=comm4d.capture(nfr)
    print('Capture completed!')
    print('Remember to transfer the webDAQ file')
    print(tn)
    time.sleep(3)
    rp.set_wave1(0.01, 0, 10, 'SINE')
    
    #os.system("read -p 'Press Enter to produce the 4D frames...' var")
    sp.wdsync()
    wdfile=sp.lastwdfile()
    print(wdfile)
    comm4d.produce(tn)
    comm4d.frames_transfer(tn)
    comm4d.save_acqdata(tn,freq4d, ampPI, freqPI)
    sp.freq4d=freq4d
    print('Freq 4d:')
    print(sp.freq4d)
    if bound is not None:
        boundfreq=bound
    else:
        boundfreq=np.array([freqPI-2, freqPI+2])

    #os.system("read -p 'Press Enter to run the analysis...' var")
    #wdfile=sp.lastwdfile()
    #print(tn)
    #print(wdfile)
    if process is None0:
        slog.log_data([tn, wdfile],[freq4d,ampPI])
    else:
        sp.peak_analysis(tn, wdfile=wdfile,bound=boundfreq)
    print(tn)
    print(wdfile)

def acq_pulse(ampPI, freqPI):
    freq4d  =264
    #rp.waves_on()
    #time.sleep(1)
    #rp.trigg4d_pulse(freq4d)
    #time.sleep(0.2)
    #rp.pulse_train(1,0.1,0,50,0.3,2)
    #rp.pulse_train(1,ampPI,0,freqPI,0.1)
    nfr = nframes
    rp.splatt_pulse(freqPI,freq4d,ampPI)
    #rp.splatt_trigger(freqPI, freq4d, ampPI)
    #rp.splatt_trigger2(freqPI, freq4d, ampPI)  # with pulse
    print('Command to RedPitaya: done!')
    print('Now waiting for transitory oscillations to damp')
    time.sleep(5)
    webdaq.start_schedule()

    tn=comm4d.capture(nfr)
    print('Capture completed!')
    print('Remember to transfer the webDAQ file')
    print(tn)
    rp.set_wave1(0.01, 0, 10, 'SINE')
    #rp.clear_wave(1)

    #os.system("read -p 'Press Enter to produce the 4D frames...' var")
    comm4d.produce(tn)
    sp.wdsync()
    comm4d.frames_transfer(tn)
    comm4d.save_acqdata(tn,freq4d, ampPI, freqPI)
    '''
    sp.freq4d=freq4d
    print('Freq 4d:')
    print(sp.freq4d)
    if bound is not None:
        boundfreq=bound
    else:
        boundfreq=np.array([freqPI-2, freqPI+2])

    os.system("read -p 'Press Enter to run the analysis...' var")
    wdfile=sp.lastwdfile()
    print(tn)
    print(wdfile)
    sp.peak_analysis(tn, wdfile=wdfile,bound=boundfreq)
    '''
    wdfile=sp.lastwdfile()
    print(tn)
    print(wdfile)

def acq_bench(ampPI, freqPI):
    rp.splatt_trigger2(freqPI, 264, ampPI)  # with pulse
    pauseStart=5
    print('Command to RedPitaya: done!')
    print('Now waiting for transitory oscillations to damp')
    time.sleep(pauseStart)
    webdaq.start_schedule()
    time.sleep(10)
    sp.wdsync()
    time.sleep(1)
    rp.set_wave1(0.01, 0, 10, 'SINE')
    #os.system("read -p 'Press Enter after sync the WD file... ' var")

    wdfile = sp.lastwdfile()
    peakOBBint, peakStandint, peakAcc2int, peakAcc3int=sp.acc_peak_analysis(wdfile,bound=np.array([freqPI-2, freqPI+2]), allaccel=1)
    print(peakOBBint, peakStandint, peakAcc2int, peakAcc3int)
    slog.log_data(['NA', wdfile],[0,ampPI,freqPI, peakOBBint, peakStandint, peakAcc2int, peakAcc3int ])

def acq_bench_seq(amp):
    startFreq= 20
    endFreq = 120
    spa = 10
    npo = np.int((endFreq-startFreq)/spa+1)
    #amp = 2
    fvec = np.linspace(startFreq, endFreq, npo)
    for i in fvec:
        s = ('Now running at: %s Hz' %i)
        print(s)
        acq_bench(amp, i)

def acq_sweep(nframes, freq4d):
    
    #suggested nframes 6000
    #  
    print('Connect Tektronik to 4D and Piezo, RedPitaya#2 to Tek trigger')
    rp.single_pulse(2)
    time.sleep(0.1)
    webdaq.start_schedule()

    tn=comm4d.capture(nframes)
    print('Capture completed!')
    print(tn)
    time.sleep(1)
    sp.wdsync()
    wdfile=sp.lastwdfile()
    print(wdfile) 
    comm4d.produce(tn)
    comm4d.frames_transfer(tn)
    comm4d.save_acqdata(tn,freq4d, 0, 0)
    slog.log_data([tn, wdfile,'Sweep'],[freq4d])

    print(tn)
    print(wdfile)

def test_sweep():
    #rp.clear_wave(2)
    #time.sleep(0.2)
    #rp.wave_on(2)
    rp.single_pulse(2)
    t0=time.time()
    #rp.splatt_pulse(1,250,0.1)
    time.sleep(0.1)
    webdaq.start_schedule()
    t1=time.time()
    print(t1-t0)
    time.sleep(15)
    sp.wdsync()
    wdfile=sp.lastwdfile()
    spea, fa=sp.acc_spectrum(wdfile)
    plot(fa, spea[0,:])
    #yscale('log')
    #xlim(0,200)
