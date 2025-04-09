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

# Define piston/tip/tilt
coords = dm.actCoords
x,y = coords[:,0], coords[:,1]
piston = np.ones(dm.nActs)
tip = x/np.max(np.abs(coords))
tilt = y/np.max(np.abs(coords))

# Define useful variables
amp = 10e-6
delay = 0.5

cmds = np.vstack((tip,tilt,piston)) 
cmds = np.array(cmds)

# Define frequency range
freq_vec = np.arange(10,130,12)

wdf_list = []
tn_list = []

for cmd in cmds:

    cmd = cmd.T * amp
    cmd = cmd.tolist()

    for k, freq in enumerate(freq_vec):

        print(f'Testing frequency: {freq:1.2f} [Hz]')
        buf_tn = eng.read(f'prepareDynCmdHistory({cmd},{freq})')

        webdaq.start_schedule()
        eng.send(f'pause({delay}); sendCmdHistory(buffer)')
        webdaq.stop_schedule()

        sp.wdsync()
        wdfile = sp.last_wdfile()
        data = sp.openfile(wdfile)
        sp.plot_data(data,ch_ids = np.array([0,1],dtype=int))

        wdf_list.append(wdfile)
        tn_list.append(buf_tn)

    print(wdf_list)
    print(tn_list)



