# import numpy as np
# import matplotlib.pyplot as plt

import time
from datetime import datetime

from devices.webDAQ import WebDAQ as wbdq
# from devices.moxa_io import Moxa_ai0
import acceleration_analysis as sp
from devices.wavegenerators import WaveGenerator

# import opticalib
# pyconf = '/mnt/libero/SPLATTData/Data/SysConfigurations/configuration.yaml'
# opticalib.load_configuration_file(pyconf)

# Connect to WebDAQ
webdaq = wbdq()
webdaq.connect()

# Start accelerometers
print('Starting accelerometers to remove startup transient')
webdaq.start_schedule()
time.sleep(10)
webdaq.stop_schedule()

# Connect to wavegenerator
wg = WaveGenerator()

def perform_sweep(fmin,fmax,ampCh1=1,ampCh2=1,duration=10):

    wg.set_wave(ch=1,ampl=ampCh1,offs=0,freq=fmin,wave_form='SIN')
    wg.set_wave(ch=2,ampl=ampCh2,offs=0,freq=fmin,wave_form='SIN')
    time.sleep(4) # wait for steady state
    wg.set_sweep(1,fmin,fmax,duration,amp=ampCh1)
    wg.set_sweep(2,fmin,fmax,duration,amp=ampCh2)

    # Start acquisition and sweep
    webdaq.start_schedule()
    wg.trigger_sweep(chPI)

    # Switch sweep mode off
    wg.sweep_mode_off(chPI)
    wg.waves_off()

    webdaq.stop_schedule_when_job_ends(job_id = 0) # Wait for webDAQ acquisition

    tnnow = datetime.now()
    tn = tnnow.strftime('%Y%m%d_%H%M%S')
    sp.wdsync()
    wdf = sp.last_wdfile()
    sp.savefile(wdf, tn)
    return tn

tn = perform_sweep(fmin=1500,fmax=8000)
