import numpy as np
import matplotlib.pyplot as plt

import time

from splattsw.devices.webDAQ import WebDAQ as wbdq
from splattsw.devices.powersupplier import PowerSupplier
from splattsw import acceleration_analysis as sp

import Pyro4
# from splattsw.devices.matlab_engine import MatlabEngine as eng

import os


# Connect to the engine
eng = Pyro4.Proxy("PYRO:matlab_engine@193.206.155.220:9090")
eng.start_engine()

# Connect to WebDAQ
webdaq = wbdq()

# Connect to power supplier
ps = PowerSupplier()
ps.load_default_state()
ps.switch_on(1)
ps.switch_on(2)
ps.switch_on(3)

# Perform  init
eng.send_command('splattInit')
eng.send_command('splattStartup')

# Set the shell
eng.send_command('splattFastSet(100e-6)')

eng.send_command('modalBase = sys_data.ff_v;')
eng.send_command('force_amp = max(2e+3,sys_data.ff_w/10);')

Nmodes = 3

# Connect to WebDAQ
webdaq.connect()

wdf_list = []
tn_list = []

for k in range(Nmodes):
    webdaq.start_schedule()

    tn = eng.get_data('splattForceStepResponse(modalBase(:,'+str(k+1)+'),force_amp('+str(k+1)+'))',n_args_out=1,is_numeric=False)
    tn_list.append(tn)

    sp.wdsync() # does os.system('rsync -av '+ftpwebdacq+' '+basepathwebdaq)
    wdfile = sp.last_wdfile()
    data = sp.openfile(wdfile)
    sp.plot_data(data,N_ch=2)

    wdf_list.append(wdfile)

print(wdf_list)

# Dock the shell
eng.send_command("splattRIP")

# Switch off all channels
ps.switch_off(3)
ps.switch_off(1)
ps.switch_off(2)

# Kill the engine
eng.stop_engine()
