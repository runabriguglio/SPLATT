from SPLATT.splattsw import splatt_analysis as sp
from SPLATT.splattsw.devices.webDAQ import WebDAQ as wdq
from SPLATT.splattsw.devices import redpitaya as rp
from matplotlib.pyplot import *

import os
import numpy as np
import time
webdaq= wdq()
webdaq.connect()

freq = 100
rp.waves_off()
rp.set_wave1(0.1,0.,freq,'SINE')
rp.set_wave2(0.1,0.,freq,'SINE')
rp.waves_on()
rp.set_phase(2,45)

fwd=1651
frp = fwd/6
bound = np.array([frp-2,frp+2])
ph=0

def acqtest(freq,ph):
    rp.waves_off()
    rp.set_wave1(0.1,0.,freq,'SINE')
    rp.set_wave2(0.1,0.,freq,'SINE')
    #rp.waves_on()
    rp.set_phase(2,ph)
    rp.waves_on() #the order of the ON sequence does not change the found phase!!

    webdaq.start_schedule()

def processtest(tfreq):
    wdfile = sp.lastwdfile()
    wds=sp.openfile(wdfile)
    clf()
    plot(wds[2,:],'-x')
    plot(wds[3,:],'-x')
    xlim(0,50)
    spe, f      = sp.acc_spectrum(wds)
    bound = np.array([tfreq-2,tfreq+2])
    peak1, peakfreq1 = sp.find_peak(spe[2,:], freq=f, bound=bound)
    peak2, peakfreq2 = sp.find_peak(spe[3,:], freq=f, bound=bound)
    phase1    = sp.peak_phase(wds[2,:], f, peakfreq1)/np.pi*180
    phase2  = sp.peak_phase(wds[3,:], f, peakfreq2)/np.pi*180
    print('Rel Phase:')
    print(phase2-phase1)
    print('Peak1')
    print(peak1)
    print('Peak2')
    print(peak2)
    print('peakfreq1')
    print(peakfreq1)
    print('peakfreq2')
    print(peakfreq2)





