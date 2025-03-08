import numpy as np
import matplotlib.pyplot as plt

import time
import os

from splattsw.devices.webDAQ import WebDAQ as wbdq
from splattsw.devices.powersupplier import PowerSupplier
from splattsw.devices.moxa_io import Moxa_ai0
import splatt_utilities as utils
from splattsw import acceleration_analysis as sp
from splattsw.devices.wavegenerators import WaveGenerator

eng = utils.start_matlab_engine()

# Connect to WebDAQ
webdaq = wbdq()
webdaq.connect()

# Start accelerometers
print('Starting accelerometers to remove startup transient')
webdaq.start_schedule()
job_status = webdaq.get_jobs_status()
while job_status[0] != 'completed':
    time.sleep(10)
    job_status = webdaq.get_jobs_status()
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
eng.send_command('splattInit')
eng.send_command('splattStartup')

# Set the shell
eng.send_command('splattFastSet()')

# Define useful data
V = eng.get_data('sys_data.ff_v')
dec = 2
eng.send_command(f'clear opts; opts.dec = {dec:%d}; opts.save2fits = 1; opts.save2mat = 0; opts.sampleNr = 256; opts.saveCmds = 1')

# Define frequency range
freq_vec = np.arange(10,130,10)

wdf_list = []
tn_list = []

for j,freq in enumerate(freq_vec):

    # Send wave
    wg.set_wave1(amp,0,freq,'SIN')
    time.sleep(2) # wait for piezo command

    # Start buffer acquisition
    eng.send_command("[pos,cur,tn]=splattAcqBufInt({'sabi32_Distance','sabi32_pidCoilOut'},opts)")

    # Start WebDAQ acquisition
    webdaq.start_schedule()

    # Wait for buffer acquisition to end
    time.sleep(30)
    tn = eng.get_data('tn',is_numeric=False)
    tn_list.append(tn)

    # Wait for webdaq acquisition to end
    job_status = webdaq.get_jobs_status()
    while job_status[0] != 'completed':
        time.sleep(1)
        job_status = webdaq.get_jobs_status()
    webdaq.stop_schedule()

    sp.wdsync()
    wdfile = sp.last_wdfile()
    # data = sp.openfile(wdfile)
    # sp.plot_data(data,ch_ids = np.array([1,3],dtype=int))

    wdf_list.append(wdfile)


print(wdf_list)
print(tn_list)
print(pres)

# Dock the shell
eng.send_command("splattRIP")

# Switch off all channels
ps.switch_off(3)
ps.switch_off(1)
ps.switch_off(2)

# Kill the engine
eng.stop_engine()

# Post processing
for k,buf_tn in enumerate(tn_list):

    data, data_addr, t_vec = utils.read_buffer_data(buf_tn)
    pos = data[0]
    dpos = pos - np.mean(pos,axis=1)
    dmode = V.T @ dpos

    modal_spe, f_vec = utils.spectral_analysis(dmode, dt = (dec+1)/1818)

    plt.figure()
    plt.plot(f_vec,modal_spe[:3])
    plt.grid('on')
    plt.legend('Mode 1','Mode 2','Mode 3')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Position [counts]')
    plt.title(f'Piezo frequency: {freq_vec[k]:%d} [Hz]')
    plt.show()

