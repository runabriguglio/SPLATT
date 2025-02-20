import numpy as np
import matplotlib.pyplot as plt

import time
import os

from splattsw.devices.webDAQ import WebDAQ as wbdq
from splattsw.devices.powersupplier import PowerSupplier
from splattsw.devices.moxa_io import Moxa_ai0
from splatt_utilities import start_matlab_engine
from splattsw import acceleration_analysis as sp

eng = start_matlab_engine()

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

# Perform  init
eng.send_command('splattInit')
eng.send_command('splattStartup')

# Set the shell
eng.send_command('splattFastSet()')

eng.send_command('modalBase = sys_data.ff_v;')
eng.send_command('force_amp = max(2e+3,sys_data.ff_w/10);')

Nmodes = 3
Nit = 6

wdf_list = []
tn_list = []

for j in range(Nit):

    if np.max((0,j-1))%2:
        eng.send_command('splattMoveBy(100e-6)')

    for k in range(Nmodes):
        webdaq.start_schedule()

        tn = eng.get_data('splattForceStepResponse(modalBase(:,'+str(k+1)+'),force_amp('+str(k+1)+'))',n_args_out=1,is_numeric=False)
        tn_list.append(tn)

        # Wait for webdaq acquisition to end
        job_status = webdaq.get_jobs_status()
        while job_status[0] != 'completed':
            time.sleep(1)
            job_status = webdaq.get_jobs_status()
        webdaq.stop_schedule()

        sp.wdsync()
        wdfile = sp.last_wdfile()
        data = sp.openfile(wdfile)
        sp.plot_data(data,ch_ids = np.array([1,3],dtype=int))

        wdf_list.append(wdfile)


# wdf_list = np.array(wdf_list).reshape([int(Nit/2),2,Nmodes])
# tn_list = np.array(tn_list).reshape([int(Nit/2),2,Nmodes])

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
