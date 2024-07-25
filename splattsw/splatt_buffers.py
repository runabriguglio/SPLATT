import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits
import os
import glob
from SPLATT.splattsw import splatt_plot as splt



SPLATT_BUFFER_FOLDER = '/home/labot/Desktop/Data/SPLATT/Buffer'

def read_buffer_data(TN = None):

    if TN is not None:
        where = os.path.join(SPLATT_BUFFER_FOLDER,TN)
    else:
        buffer_folder_list = sorted(glob.glob(SPLATT_BUFFER_FOLDER))
        where = buffer_folder_list.split("/")[-1]

    which = os.path.join(where,'Data.fits')
    hdu = pyfits.open(which)
    raw_data = hdu[0].data
    raw_data_size = np.shape(raw_data)

    if len(raw_data_size) == 3:
        time_vec = raw_data[0][0]
        pos = raw_data[0][1:]
        force = raw_data[1][1:]
        data = np.array([pos,force])
    else:
        time_vec = raw_data[0]
        data = raw_data[1:]

    return data, time_vec


def analyse_buffer_data(TN = None, show = False):

    data, time_vec = read_buffer_data(TN)

    data_size = np.shape(data)
    dt = time_vec[1]-time_vec[0]

    if len(data_size) == 3:
        max_osc = np.zeros(data_size[:-1])
        peak_freq = np.zeros(data_size[:-1])
        for k in range(data_size[0]):
            spe, f = spectral_analysis(data[k],dt)

            if show:
                plot_freq_data(f,spe)
                plot_data(time_vec,data[k])

            max_osc[k], peak_freq[k] = find_peak_freq(spe,f,bound=[1.,f[-1]])
            splt.plot_splatt_data(max_osc[k])
    else:

        spe, f = spectral_analysis(data,dt)
        spe_size = np.shape(spe)
        max_osc = np.zeros(spe_size)
        peak_freq = np.zeros(spe_size)

        if show:
            plot_freq_data(f,spe)
            plot_data(time_vec,data)

        max_osc, peak_freq = find_peak_freq(spe,f,bound=[1.,f[-1]])
        splt.plot_splatt_data(max_osc)

    return max_osc, peak_freq


def analyse_oscillation(TN_list, freq_list):

    peak_val = np.zeros([19,len(TN_list)])
    peak_freq = np.zeros([19,len(TN_list)])

    for k,TN in enumerate(TN_list):
        freq = freq_list[k]
        data, t_vec = read_buffer_data(TN)
        dt = t_vec[1] - t_vec[0]
        spe, f_vec = spectral_analysis(data,dt)

        if freq is not None:
            freq_bound = [freq - 2., freq + 2.]
        else:
            freq_bound = [1., f_vec[-1]]

        peak_v, peak_f = find_peak_freq(spe, f_vec, freq_bound)
        splt.plot_splatt_data(peak_v)

        # Store data in output variables
        peak_val[:,k] = peak_v
        peak_freq[:,k] = peak_f

        # remove mean from data
        data_m = np.mean(data,axis=1)
        data_m = np.repeat(data_m,len(data[0,:]))
        data_mean = data_m.reshape(np.shape(data))
        data_osc = data - data_mean
        data_osc = data_osc/(2.**26)
        x = splt.get_splatt_act_coord()[:,0]
        y = splt.get_splatt_act_coord()[:,1]
        x_rep = np.ones([19,1])*x
        y_rep = np.repeat(y,19,axis=0)
        x_rep = x_rep.reshape([19,19])
        y_rep = y_rep.reshape([19,19])

        x_coeffs = x_rep @ data_osc
        y_coeffs = y_rep @ data_osc

        print(x_coeffs)
        print(y_coeffs)

    return peak_val, peak_freq





def plot_data(tvec, data):

    for ind in range(data[0]):
        plt.figure()
        plt.scatter(tvec,data)
        plt.xlabel('Time [s]')
        plt.ylabel('Amplitude [bits]')
        plt.title('Coil ' + str(ind))
        plt.show()


def plot_freq_data(freq, data):

    for ind in range(data[0]):
        plt.figure()
        plt.scatter(freq,data)
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [bits]')
        plt.title('Coil ' + str(ind))
        plt.show()


def spectral_analysis(signal, dt = 1):
    """ This is equivalent to:
        from M4/m4/mini_OTT import time_history as th
        th.spectrum(signal,dt)
        """

    print(dt)
    nsig = signal.shape
    if np.size(nsig) == 1:
        thedim = 0
    else:
        thedim = 1

    if thedim == 0:
        spe = np.fft.rfft(signal, norm="ortho")
        nn = np.sqrt(spe.shape[thedim])
    else:
        spe = np.fft.rfft(signal, axis=1, norm="ortho")
        nn = np.sqrt(spe.shape[thedim])

    spe = (np.abs(spe)) / nn
    freq = np.fft.rfftfreq(signal.shape[thedim], d=dt)

    if thedim == 0:
        spe[0] = 0
    else:
        spe[:, 0] = 0

    return spe, freq


def find_peak_freq(spe, freq_vec, bound = None):
    """ This differs from the find_peak() in SPLATT/splattsw/splatt_analysis.py
        as it normalizes the amplitude to the signal amplitude and has the frequency
        vector as mandatory input. As a result, it also works for 2D arrays
        """
    freq_len = len(freq_vec)
    idf = np.arange(freq_len)

    if bound is not None:
        idf = np.where(np.logical_and(freq_vec>= bound[0], freq_vec<=bound[1]))

    size_spe = np.shape(spe)
    nChan = size_spe[0]

    if freq_len is not size_spe[1]:
        spe = spe.T
        nChan = size_spe[0]

    peak_val = np.zeros(nChan)
    peak_freq = np.zeros(nChan)

    for j in range(nChan):
        v = spe[j]
        peak = np.max(v[idf])
        peak_val[j] = peak*np.sqrt(2) # re-normalize amplitude
        peak_id = np.argmax(v[idf])
        idf1 = idf[0]
        peak_freq[j] = freq_vec[idf1[peak_id]]

    return peak_val, peak_freq