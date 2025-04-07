import numpy as np
import matplotlib.pyplot as plt

import time
import os

from devices.webDAQ import WebDAQ as wbdq
from devices.powersupplier import PowerSupplier
#from devices.moxa_io import Moxa_ai0
import splatt_utilities as utils
import acceleration_analysis as sp
from devices.wavegenerators import WaveGenerator
from devices.deformable_mirror import SPLATTEngine

# Connect to WebDAQ
webdaq = wbdq()
webdaq.connect()

# Start accelerometers
print('Starting accelerometers to remove startup transient')
webdaq.start_schedule()
time.sleep(10)
webdaq.stop_schedule()

dm = SPLATTEngine()
eng = dm._eng

# # Connect to moxa
# mx = Moxa_ai0()
# pres = mx.read_pressure()

# # Connect to power supplier
# ps = PowerSupplier()
# ps.load_default_state()
# time.sleep(1) # wait for load
# ps.switch_on(1)
# ps.switch_on(2)
# ps.switch_on(3)

# # Connect to wavegenerator
# wg = WaveGenerator()
# amp = 1

# # Perform  init
# eng.send('splattStartup')

# # Set the shell
# eng.send('splattFastSet()')

# # Define useful data
# V = eng.read('sys_data.ff_v')
# dec = 2
# eng.send(f'clear opts; opts.dec = {dec:%d}; opts.save2fits = 1; opts.save2mat = 0; opts.sampleNr = 256; opts.saveCmds = 1')

# # Define frequency range
# freq_vec = np.arange(10,130,10)

wdf_list = []
# tn_list = []

Nrep = 3
amp_vec = np.array([1000,1500,2000])

for kk in range(Nrep):
    fcmd = np.array(eng.read("aoRead('sabi16_force',1:19)"))
    fcmd = np.reshape(fcmd,dm.nActs)
    old_fcmd = fcmd.tolist()

    # Start WebDAQ acquisition
    webdaq.start_schedule()

    time.sleep(2) # wait for schedule to actually start

    for j,amp in enumerate(amp_vec):
        new_fcmd = (old_fcmd+amp).tolist()
        eng.send(f'lattApplyForce({new_fcmd})')
        eng.send(f'lattApplyForce({old_fcmd})')
        print(f'{j+1}/{len(amp_vec)}')

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

    wdf_list.append(wdfile)

print(wdf_list)


# Dock the shell
eng.send("splattRIP")

# Switch off all channels
ps.switch_off(3)
ps.switch_off(1)
ps.switch_off(2)

# Kill the engine
eng.stop_engine()

# Load and normalize
V = np.loadtxt('/home/labot/git/SPLATT/SPLATT_Data/mirror_modes.txt')
V = V/np.sqrt(np.shape(V)[0])

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

