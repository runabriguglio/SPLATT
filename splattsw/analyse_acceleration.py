import splatt_analysis as sp

import matplotlib.pyplot as plt
import numpy as np

sp.wdsync()
# wdfile = 'OBB-Vibration_2024-12-09T1036-55-939.wdd'
# wdfile='OBB-Vibration_2024-12-03T17-28-03-210.wdd' # during step, no FF
wdfile = 'OBB-Vibration_2024-12-05T11-25-08-886.wdd' # during step, with FF

n_samples = int(3e+6) # number of samples (for manual acquisition only)
n_ch = 2 # number of channels

data = sp.openfile(wdfile, data_len = n_samples)
sp.plot_data(data, N_ch = n_ch)

data_ch0 = data[0]
spe0, f0 = sp.acc_spectrum(data_ch0)
peak, peak_freq0, peak_id = sp.find_peak(spe0,f0)
amp0 = sp.acc_integrate(spe0, peak_freq0, peak_id)

bound = np.array([peak_freq0-10,peak_freq0+10])
data_ch1 = data[1]
spe1, f1 = sp.acc_spectrum(data_ch1)
peak, peak_freq1, peak_id = sp.find_peak(spe1,f1,freq_bounds=bound)
amp1 = sp.acc_integrate(spe1, peak_freq1, peak_id)

# net acceleration
net_acc = data_ch0 - data_ch1
spe, f = sp.acc_spectrum(net_acc)
peak, peak_freq, peak_id = sp.find_peak(spe,f)
amp = sp.acc_integrate(spe, peak_freq, peak_id)

plt.figure()
plt.plot(net_acc)
plt.grid('on')
plt.ylabel('Amplitude [g]')
plt.title('Acceleration difference')
plt.show()

plt.figure()
plt.plot(f,spe)
plt.grid('on')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [g]')
plt.title('Acceleration difference spectrum')
plt.show()
