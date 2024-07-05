from SPLATT.splattsw import splatt_analysis as sp
from SPLATT.splattsw.devices.webDAQ import WebDAQ as wbdq
from SPLATT.splattsw.devices import devices as wg

import matlotlib.pyplot as plt
import numpy as np
from time import sleep

# Connect to WebDAQ
webdaq = wbdq()
webdaq.connect()

freqV = np.arange(510,120,5)
N = len(freqV)
Nch = 4
maxD = np.zeros([Nch,N],dtype=float)
maxF = np.zeros([Nch,N],dtype=float)

ampPI = 1
ampGain = 1

# Select device ('Rigol' or 'RedPitaya')
wg.update_device('Rigol')

for i in range(N):

    freqPI = freqV[i]

    wg.set_wave1(ampPI/ampGain,0,freqPI,'SINE')
    time.sleep(2) # wait for piezo command
    webdaq.start_schedule()
    time.sleep(4) # wait for acquisition to end

    #accelerometer analysis
    #read file
    sp.wdsync() # does os.system('rsync -av '+ftpwebdacq+' '+basepathwebdaq)
    wdfile = sp.lastwdfile('Piezo')

    #spectrum
    spe, f = sp.acc_spectrum(wdfile)
    # sp.accplot(f, spe[0,:])

    # Loop on all channels
    for j in range(Nch):
        peak, maxF[j,i]  = sp.find_peak(spe[j,:],f)
        maxD[j,i] = sp.acc_integrate(peak,maxF[j,i])

#print(peaks)
for chId in range(Nch):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(maxF[chId,:],maxD[chId,:])
    ax.plot(maxF[chId,:],maxD[chId,:],"o")
    ax.set_xlabel("Peak frequency [Hz]")
    ax.set_ylabel("Peak oscillation [um]")
    #ax.set_yscale("log")



