from SPLATT.splattsw import splatt_analysis as sp
from SPLATT.splattsw.devices.webDAQ import WebDAQ as wbdq
from SPLATT.splattsw.devices import devices as wg

import matplotlib.pyplot as plt
import numpy as np
import time


# Connect to WebDAQ
webdaq = wbdq()
webdaq.connect()
freq_vec = np.arange(20,100+5,5)


def analyse_latest_wdfile(doplot = False, chId = 0):
    #accelerometer analysis
    #read file
    sp.wdsync() # does os.system('rsync -av '+ftpwebdacq+' '+basepathwebdaq)
    wdfile = sp.lastwdfile('Piezo')

    #spectrum
    spe_4ch, f = sp.acc_spectrum(wdfile)
    spe = spe_4ch[chId,:]

    if doplot == True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(f,spe)
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel("Amplitude [g]")
        fig.show()

    peak, peak_f = sp.find_peak(spe,freq=f,bound=[1.,f[-1]])
    peak_d = sp.acc_integrate(peak,peak_f)

    return peak_d, peak_f

def test_single_freq(freq_val,amp=1):
    wg.set_wave1(amp,0,freq_val,'SINE')
    time.sleep(2) # wait for piezo command
    webdaq.start_schedule()
    time.sleep(40) # wait for acquisition to end
    print('Done!')


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
        test_single_freq(freqPI) #ampPi/ampGain
        maxD[i], maxF[i] = analyse_latest_wdfile()

    plt.figure()
    plt.plot(maxF, maxD)
    plt.scatter(maxF, maxD, marker='o', c='red', s=15)
    plt.xlabel("Peak frequency [Hz]")
    plt.ylabel("Peak oscillation [m]")
    plt.show()

    return maxD, maxF



