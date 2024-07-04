#memo to run the SPLATT measurements in the laboratory
#importing
from SPLATT.splattsw import splatt_analysis as sp
#from SPLATT.splattsw import splatt_acq as acq

#from importlib  import reload #for reload
#power on and connect the RedPitaya
#start the redpitaya server
#systemctl start redpitaya_scpi
from SPLATT.splattsw.devices.webDAQ import WebDAQ as wbdq
from SPLATT.splattsw.devices import redpitaya as rp

import numpy as np
import time
#import os

webdaq = wbdq()
webdaq.connect()

# end of initialization commands
#see cmdlog.py for the commands to be used

freqV = np.linspace(10,120,20)
peaks = np.zeros([len(freqV),1],dtype=float)

freqPI = 10
ampPI = 2
ampGain = 10

for i in range(len(freqV)):

   webdaq.start_schedule()

    # INSERT TIME SLEEP?

   rp.set_wave1(ampPi/ampGain,0,freqPi,'SQUARE')
   time.sleep(6)
   rp.clear_wave1()
   webdaq.stop_schedule()
   os.system('rsync -av '+ftpwebdacq+' '+basepathwebdaq)
   wdfile = sp.lastwdfile()

    #accelerometer analysis
    #read file
    #spectrum
    spe, f = sp.acc_spectrum(wdfile)
    p.accplot(f, spe[0,:])

    peak, peakf = sp.find_peak(spe,f)
    maxd = sp.acc_integrate(peak,peakf)

   peaks[i] = maxd

print(peaks)
