from SPLATT.splattsw import splatt_analysis as sp
from SPLATT.splattsw.devices.webDAQ import WebDAQ as wbdq
from SPLATT.splattsw.devices import devices as wg

import matplotlib.pyplot as plt
import numpy as np
import time


# Connect to WebDAQ
webdaq = wbdq()
webdaq.connect()
freq_vec = np.arange(40,50+1,1)

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

    return peak, peakfreq, idf1[peakid]

def acc_integrate(spe, peak_freq, peak_id, delta_peak = 3):
    #acc is the acceleration value in g
    #print('Is the acceleration provided as [g]?')
    spe_peak = spe[(peak_id-delta_peak):(peak_id+delta_peak+1)]
    acc = np.sum(spe_peak)
    acc = acc * 9.807#converts to m/s2
    amp = acc/(4*np.pi**2*peak_freq**2)

    return amp


def analyse_wdfile(wdfile, doplot = False, ch = 0):
    #accelerometer analysis
    #spectrum
    spe_4ch, f = sp.acc_spectrum(wdfile)
    spe = spe_4ch[ch,:]*np.sqrt(2) # re-multiply by sqrt(2) to have amp=oscillation_amp

    if doplot == True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(f,spe)
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel("Amplitude [g]")
        fig.show()

    peak, peak_f, peakId = find_peak(spe,freq=f,bound=[1.,f[-1]])#sp.find_peak(spe,freq=f,bound=[1.,f[-1]])
    peak_d = acc_integrate(spe,peak_f,peakId) #sp.acc_integrate(peak,peak_f)

    return peak_d, peak_f


def test_single_freq(freq_val,amp=1):
    wg.set_wave1(amp,0,freq_val,'SINE')
    time.sleep(2) # wait for piezo command
    webdaq.start_schedule()
    time.sleep(40) # wait for acquisition to end
    print('Done!')

    #read file
    sp.wdsync() # does os.system('rsync -av '+ftpwebdacq+' '+basepathwebdaq)
    wdfile = sp.lastwdfile('Piezo')

    return wdfile


def start_cycle_on_freq(freqV = freq_vec):
    ampPI = 1
    ampGain = 1

    # Select device ('Rigol' or 'RedPitaya')
    wg.update_device('Rigol')

    N = len(freqV)
    maxD = np.zeros(N,dtype=float)
    maxF = np.zeros(N,dtype=float)

    for i in range(N):
        freqPI = freqV[i]
        wdf=test_single_freq(freqPI) #ampPi/ampGain
        maxD[i], maxF[i] = analyse_wdfile(wdf)

    plt.figure()
    plt.plot(maxF, maxD)
    plt.scatter(maxF, maxD, marker='o', c='red', s=15)
    plt.xlabel("Peak frequency [Hz]")
    plt.ylabel("Peak oscillation [m]")
    plt.show()

    return maxD, maxF









