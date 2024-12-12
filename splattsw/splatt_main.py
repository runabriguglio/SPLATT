import numpy as np
import matplotlib.pyplot as plt

import Pyro4
import time

from devices.webDAQ import WebDAQ as wbdq
from devices import powersupplier as ps

# Connect to WebDAQ
webdaq = wbdq()

# Connect to the engine
eng = Pyro4.Proxy("PYRO:matlab_engine@193.206.155.220:9090")
eng.start_engine()

# Connect to power supplier
ps.load_saved_state()
ps.switch_on(1)
ps.switch_on(2)

# Perform  init
eng.send_command('splattInit')

# Test commands
for k in range(5):
    str_text = 'disp(' + str(k*10+1) + ')'
    eng.send_command(str_text)

# # Switch on coils
# ps.switch_on(3)
# eng.send_command('splattStartup')
#
# # Connect to WebDAQ
# webdaq.connect()
#
# webdaq.start_schedule()
# time.sleep(time_for_acquisition) # wait for acquisition to end
#
#
# # Kill the engine
# eng.stop_engine()
