from SPLATT.splattsw import splatt_analysis as sp

import matplotlib.pyplot as plt
import numpy as np

sp.wdsync()
wdfile = sp.last_wdfile()

n_samples = int(3e+5) # number of samples (for manual acquisition only)
n_ch = 2 # number of channels

data = sp.openfile(wdfile, data_len = n_samples)
sp.plot_data(data, N_ch = n_ch)

data_ch0 = data[0]
spe0, f0 = sp.acc_spectrum(data_ch0)
peak, peak_freq, peak_id = sp.find_peak(spe0)
amp = sp.acc_integrate(spe0, peak_freq, peak_id)