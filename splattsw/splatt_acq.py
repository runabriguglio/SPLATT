
#from splattsw import splatt_analysis as sp
from splattsw import acceleration_analysis as sp
from splattsw.devices.wavegenerators import WaveGenerator as wg
from splattsw.devices.webDAQ import WebDAQ as wdq

#import aoptics
#pyconf = '/mnt/libero/SPLATTData/Data/SysConfigurations/configuration.yaml'
#aoptics.load_configuration_file(pyconf)

from aoptics.devices.interferometer import PhaseCam4020
from splattsw.devices.moxa_io import Moxa_ai0
from splattsw.devices.deformable_mirror import SPLATTEngine

import matplotlib.pyplot as plt
import os
import time
import numpy as np

# Define devices
interf=PhaseCam4020()
webdaq= wdq()
wavegen = wg()
mx = Moxa_ai0()
eng = None

def acq_sweep(fmin = 30, fmax = 110, duration = 11, ampPI = 2, nframes:int = 2250, chPI:int = 1, produce:bool = False):


    webdaq.connect()

    # Setup sweep parameters
    wavegen.set_wave(ch=chPI,ampl=ampPI,offs=0,freq=fmin,wave_form='SIN')
    time.sleep(4)
    wavegen.sweep(chPI,fmin,fmax,duration,amp=ampPI)

    # Start acquisition and sweep
    webdaq.start_schedule()
    wavegen.trigger_sweep(chPI)
    time.sleep(0.5)
    tn=interf.capture(nframes)
    
    # Acquisition end
    #wavegen.clear_wave(chPI)
    print(f'Capture completed! Saved in: {tn}')
    webdaq.stop_schedule()

    # Post-processing
    time.sleep(1)
    sp.wdsync()
    wdfile=sp.last_wdfile()
    print(f'Acceleration data: {wdfile}')
    pres = mx.read_pressure()
    print(f'The vacuarium pressure is {pres:1.3f} [bar]')

    if produce:
        interf.produce(tn)



def acq_freq(freqPI,  ampPI=2, nframes:int = 2000, chPI:int = 1, produce:bool = False, buffer:bool = False):

    if buffer:
        # Connect to SPLATT
        
        dm = SPLATTEngine()
        eng = dm._eng
        eng.send('clear opts; opts.dec = 0; opts.sampleNR = 256; opts.save2mat = 0; opts.save2fits = 1')
    
    # Start acquisition and sine wave
    webdaq.connect()
    wavegen.set_wave(ch=chPI, ampl=ampPI, offs=0, freq=freqPI, wave_form='SIN')
    if buffer:
        eng.oneway_send("[pos,cur,tn]=splattAcqBufInt({'sabi32_Distance','sabi32_pidCoilOut'},opts)")
    t0 = time.time()
    time.sleep(4)
    webdaq.start_schedule()
    t0a = time.time()
    if nframes > 0:
        tn=interf.capture(nframes)
    
    t1a = time.time()
    dt = t1a-t0a
    if dt < 10:
        time.sleep(10-dt)

    # Acquisition end
    wavegen.wave_off(chPI)

    if nframes > 0:
        print(f'Capture completed! Saved in: {tn}')
    webdaq.stop_schedule()
    
    # Post-processing
    time.sleep(1)
    sp.wdsync()
    wdfile=sp.last_wdfile()
    print(f'Acceleration data: {wdfile}')
    p_bar = mx.read_pressure()
    print(f'The vacuarium pressure is {p_bar:1.3f} [bar]')

    if np.logical_and( nframes > 0, produce):
        interf.produce(tn)

    if buffer:
        # Wait for buffer to end before reading tn
        t1 = time.time()
        dt = t1-t0

        if dt < 45:
            time.sleep(45-dt)

        buf_tn = eng.read('tn')
        print(f'Buffer TN is: {buf_tn}')




