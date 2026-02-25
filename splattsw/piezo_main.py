# import numpy as np
# import matplotlib.pyplot as plt

import time
from datetime import datetime

from devices.webDAQ import WebDAQ as wbdq
# from devices.moxa_io import Moxa_ai0
import acceleration_analysis as sp
from devices.wavegenerators import WaveGenerator

# Connect to WebDAQ
webdaq = wbdq()
webdaq.connect()

# Connect to wavegenerator
wg = WaveGenerator()

# Start accelerometers
def start_webdaq():
    print('Starting accelerometers to remove startup transient')
    webdaq.start_schedule()
    time.sleep(10)
    webdaq.stop_schedule()

def perform_sweep(fmin,fmax,ch=1,ampCh=1,duration=10):

    wg.set_wave(ch,ampl=ampCh,offs=0,freq=fmin,wave_form='SIN')
    wg.set_sweep(ch,fmin,fmax,duration,amp=ampCh,return_time=0)

    # Start acquisition and sweep
    webdaq.start_schedule()
    time.sleep(2)
    wg.trigger_sweep(ch)

    webdaq.stop_schedule_when_job_ends(job_id = 0)

    tnnow = datetime.now()
    tn = tnnow.strftime('%Y%m%d_%H%M%S')
    sp.wdsync()
    wdf = sp.last_wdfile(ext='PI')
    sp.savefile(wdf, tn)
    return tn

#tn = perform_sweep(ch=1,fmin=1500,fmax=8000)
