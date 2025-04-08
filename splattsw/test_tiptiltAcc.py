import numpy as np
import matplotlib.pyplot as plt

import time

from devices.webDAQ import WebDAQ as wbdq
import acceleration_analysis as sp
from devices.deformable_mirror import SPLATTEngine

# Connect to SPLATT
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

# Define piston/tip/tilt
coords = dm.actCoords
x,y = coords[:,0], coords[:,1]
piston = np.ones(dm.nActs)
tip = x/np.max(np.abs(coords))
tilt = y/np.max(np.abs(coords))

# Define useful variables
amp = 3e-6
delay = 2

cmds = np.hstack((tip,tilt,piston)) * amp

# Define frequency range
freq_vec = np.arange(10,130,12)

wdf_list = []
tn_list = []

for cmd in cmds:

    cmd = cmd.tolist()

    for k, freq in enumerate(freq_vec):

        print(f'Testing frequency: {freq:1.2f} [Hz]')
        buf_tn = eng.read(f'prepareDynCmdHistory({cmd},{freq})')
        eng.oneway_send(f'pause({delay}); sendCmdHistory(buffer)')

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

        wdf_list.append(wdfile)
        tn_list.append(buf_tn)

    print(wdf_list)
    print(tn_list)



