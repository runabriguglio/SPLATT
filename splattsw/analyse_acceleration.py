import SPLATT.splattsw.utils.folder_paths
import acceleration_analysis as sp

import matplotlib.pyplot as plt
import numpy as np

sp.wdsync()
# wdfile = 'OBB-Vibration_2024-12-09T1036-55-939.wdd'
# wdfile='OBB-Vibration_2024-12-03T17-28-03-210.wdd' # during step, no FF
wdfile = 'OBB-Vibration_2024-12-05T11-25-08-886.wdd' # during step, with FF

def find_net_acceleration(wdf, n_samples = None):
    n_ch = 2 # number of channels

    data = sp.openfile(wdfile, data_len = n_samples)
    sp.plot_data(data, N_ch = n_ch)

    # data_ch0 = data[0]
    # spe0, f0 = sp.acc_spectrum(data_ch0)
    # peak, peak_freq0, peak_id = sp.find_peak(spe0,f0)
    # amp0 = sp.acc_integrate(spe0, peak_freq0, peak_id)
    #
    # bound = np.array([peak_freq0-10,peak_freq0+10])
    # data_ch1 = data[1]
    # spe1, f1 = sp.acc_spectrum(data_ch1)
    # peak, peak_freq1, peak_id = sp.find_peak(spe1,f1,freq_bounds=bound)
    # amp1 = sp.acc_integrate(spe1, peak_freq1, peak_id)

    # net acceleration
    interval=300:600
    net_acc = data[0,interval] - data[1,interval]
    spe, f = sp.acc_spectrum(net_acc)
    peak, peak_freq, peak_id = sp.find_peak(spe,f)
    amp = sp.acc_integrate(spe, peak_freq, peak_id)

    sp.plot_data(net_acc,N_ch=1,title_str='Acceleration difference: ch0-ch1')

    plt.figure()
    plt.plot(f,spe)
    plt.grid('on')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Power [g]')
    plt.title('Acceleration difference spectrum')
    plt.show()

wdlist = ['OBB-Vibration_2024-12-12T13-01-00-896.wdd',
'OBB-Vibration_2024-12-12T13-01-21-039.wdd',
'OBB-Vibration_2024-12-12T13-01-41-324.wdd',
'OBB-Vibration_2024-12-12T13-02-01-472.wdd',
'OBB-Vibration_2024-12-12T13-02-21-666.wdd',
'OBB-Vibration_2024-12-12T13-02-41-710.wdd']

for k,wdf in enumerate(wdlist):
    find_net_acceleration(wdf)


def analyse_file_list(wdlist):

    maxima = np.zeros([len(wdlist),2])
    rms_vals = np.zeros([len(wdlist),2])

    for k,file_name in enumerate(wdlist):
        data = sp.openfile(file_name)
        sp.plot_data(data,N_ch = 2)
        net_acc = data[0] - data[1]
        sp.plot_data(net_acc, N_ch=1, title_str='Acceleration difference: ch0-ch1')
        max_ch0 = np.max(np.abs(data[0]))
        max_ch1 = np.max(np.abs(data[1]))
        maxima[k,0] = max_ch0
        maxima[k,1] = max_ch1
        N = 6
        large_N0 = np.sort(np.abs(data[0]))[-N:] # 6 largest values
        large_N1 = np.sort(np.abs(data[1]))[-N:] # 6 largest values
        rms_vals[k,0] = np.sqrt(np.sum(large_N0**2)/N)
        rms_vals[k,1] = np.sqrt(np.sum(large_N1**2)/N)

    return maxima, rms_vals