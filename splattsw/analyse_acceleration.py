import splatt_analysis as sp

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
peak, peak_freq0, peak_id = sp.find_peak(spe0,f0)
amp0 = sp.acc_integrate(spe0, peak_freq0, peak_id)

data_ch1 = data[1]
spe1, f1 = sp.acc_spectrum(data_ch1)
peak, peak_freq1, peak_id = sp.find_peak(spe1,f1,bound = np.array([peak_freq0-10,peak_freq0+10]))
amp1 = sp.acc_integrate(spe1, peak_freq1, peak_id)