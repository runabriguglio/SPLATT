import numpy as np
import matplotlib.pyplot as plt

import Pyro4
import time

from splattsw.devices.webDAQ import WebDAQ as wbdq
from splattsw.devices.powersupplier import PowerSupplier
from splattsw import acceleration_analysis as sp

import os

os.system(f"ssh labot@193.206.155.220 'python3 /home/labot/git/SPLATT/splattsw/devices/matlabengine_server.py' &")

# import subprocess
#
# HOST = 'labot@193.206.155.220'
# COMMAND = 'python3 /home/labot/git/SPLATT/splattsw/devices/matlabengine_server.py'
# ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
#                            shell=True,
#                            stdout=subprocess.PIPE,
#                            stderr=subprocess.PIPE)

# Connect to the engine
eng = Pyro4.Proxy("PYRO:matlab_engine@193.206.155.220:9090")
eng.start_engine()

# Connect to WebDAQ
webdaq = wbdq()

# Connect to power supplier
ps = PowerSupplier()
ps.load_saved_state()
ps.switch_on(1)
ps.switch_on(2)

# Perform  init
eng.send_command('splattInit')
eng.send_command('splattStartup')

# Switch on coils
ps.switch_on(3)

# Set the shell
eng.send_command('splattFastSet(100e-6)')

# Connect to WebDAQ
webdaq.connect()

bias_vec = (np.arange(6)-1)*2000

wdf_list = []

for k,bias_f in enumerate(bias_vec):
    webdaq.start_schedule()

    eng.send_command("lattApplyForce("+str(bias_f)+")")
    eng.send_command("lattApplyForce(0)")
    eng.send_command("lattApplyForce("+str(bias_f)+")")
    eng.send_command("lattApplyForce(0)")
    eng.send_command("lattApplyForce("+str(bias_f)+")")
    eng.send_command("lattApplyForce(0)")

    time.sleep(5) # wait for acquisition to end

    sp.wdsync() # does os.system('rsync -av '+ftpwebdacq+' '+basepathwebdaq)
    wdfile = sp.last_wdfile()
    data = sp.openfile(wdfile)
    sp.plot_data(data,N_ch=1)

    wdf_list.append(wdfile)

print(wdf_list)

# Dock the shell
eng.send_command("splattDock")

# Switch off all channels
ps.switch_off(3)
ps.switch_off(1)
ps.switch_off(2)

# Kill the engine
eng.stop_engine()
