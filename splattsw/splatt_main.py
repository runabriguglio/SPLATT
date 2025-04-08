import numpy as np
import matplotlib.pyplot as plt

import time

from devices.webDAQ import WebDAQ as wbdq
from devices.powersupplier import PowerSupplier
from devices.moxa_io import Moxa_ai0
import acceleration_analysis as sp
from devices.wavegenerators import WaveGenerator
from devices.deformable_mirror import SPLATTEngine
import splatt_utilities as utils

dm = SPLATTEngine()
eng = dm._eng

# Connect to WebDAQ
webdaq = wbdq()
webdaq.connect()

# Start accelerometers
print('Starting accelerometers to remove startup transient')
webdaq.start_schedule()
time.sleep(10)
webdaq.stop_schedule()

# Connect to moxa
mx = Moxa_ai0()
pres = mx.read_pressure()

# Connect to power supplier
ps = PowerSupplier()
ps.load_default_state()
time.sleep(1) # wait for load
ps.switch_on(1)
ps.switch_on(2)
ps.switch_on(3)

# Connect to wavegenerator
wg = WaveGenerator()
amp = 1

# Perform  init
eng.send('splattStartup')

# Set the shell
eng.send('splattFastSet()')

# Define useful data
V = dm.mirrorModes
dec = 2
eng.send(f'clear opts; opts.dec = {dec:%d}; opts.save2fits = 1; opts.save2mat = 0; opts.sampleNr = 256')

# Define frequency range
freq_vec = np.arange(10,130,10)

wdf_list = []
tn_list = []

for k, freq in enumerate(freq_vec):

    wg.set_wave1(amp, freq=freq, wave_form='SIN')
    time.sleep(1)

    # Start buffer acquisition
    eng.oneway_send("[pos,cur,tn]=splattAcqBufInt({'sabi32_Distance','sabi32_pidCoilOut'},opts)")

    # Start WebDAQ acquisition
    webdaq.start_schedule()

    # Wait for webdaq acquisition to end
    job_status = webdaq.get_jobs_status()
    while job_status[0] != 'completed':
        time.sleep(1)
        job_status = webdaq.get_jobs_status()
    webdaq.stop_schedule()

    sp.wdsync()
    wdfile = sp.last_wdfile()
    data = sp.openfile(wdfile,18*1652)
    sp.plot_data(data,ch_ids = np.array([0,1],dtype=int))

    buf_tn = eng.read('tn')

    wdf_list.append(wdfile)
    tn_list.append(buf_tn)

print(wdf_list)
print(tn_list)

# Dock the shell
eng.send("splattRIP")

# Switch off all channels
ps.switch_off(3)
ps.switch_off(1)
ps.switch_off(2)

# Kill the engine
eng.close()

# Load and normalize
V = V/np.sqrt(np.shape(V)[0])

utils.buffsync()

# Post processing
for k,buf_tn in enumerate(tn_list):

    data, data_addr, t_vec = utils.read_buffer_data(buf_tn)
    pos = data[0]
    dpos = pos - np.mean(pos,axis=1)
    dmode = V.T @ dpos

    modal_spe, f_vec = utils.spectral_analysis(dmode, dec = dec)

    plt.figure()
    plt.plot(f_vec,modal_spe[:3])
    plt.grid('on')
    plt.legend('Mode 1','Mode 2','Mode 3')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Position [counts]')
    plt.title(f'Piezo frequency: {freq_vec[k]:%d} [Hz]')
    plt.show()

