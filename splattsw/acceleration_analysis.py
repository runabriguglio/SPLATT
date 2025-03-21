import numpy as np
from matplotlib import pyplot as plt

import os
import glob
import json
import struct

from splattsw.utils import timehistory as th

# WebDAQ variables
freqwebdaq = 1651.6129 #Hz; minimum sampling frequency
basepathwebdaq = '/mnt/jumbo/SPLATT/WebDaqData/'
#ftpwebdacq = '/home/ftpuser/ftp/' # old files
ftpwebdacq = '/home/ftpuser/ftp/files/' # new files, from 2024 Q4 onwards

def wdsync():
    os.system('rsync -av '+ftpwebdacq+' '+basepathwebdaq)

def openfile(name, data_len = None):
    file_path = os.path.join(basepathwebdaq, name)
    data = _openwdd(file_path, data_len)
    return data

def plot_data(data, ch_ids = None, freq = freqwebdaq, title_str = None):
    data_size = np.shape(data)

    if ch_ids is None:
        ch_ids = np.arange(data_size[0])

    t_vec = np.arange(0,np.max(data_size))*1/freq


    plt.figure()
    plt.legend('Channel' + str(ch_ids))
    if title_str is not None:
        plt.title(title_str)
    for ch in ch_ids:

        if len(data_size) == 1:
            vec = data
        else:
            vec = data[ch]

        plt.plot(t_vec,vec)

        plt.xlabel('Time [s]')
        plt.ylabel('Acceleration [g]')
        plt.show()
        plt.grid('on')
        plt.axis('tight')



def find_peak(v, freq=None, freq_bounds=None):

    idf = range(len(v))
    if freq is not None and freq_bounds is not None:
        idf = np.where(np.logical_and(freq>= freq_bounds[0], freq<=freq_bounds[1]))
        if len(idf) == 1:
            idf = idf[0]
    peak = max(v[idf])
    peakid = np.argmax(v[idf])

    if freq is not None:
        peakfreq = freq[idf[peakid]]
        return peak, peakfreq, idf[peakid]

    return peak, peakid

def acc_integrate(spe, peak_freq, peak_id, delta_peak = 3):
    #acc is the acceleration value in g
    spe_peak = spe[(peak_id-delta_peak):(peak_id+delta_peak+1)]
    acc = np.sum(spe_peak)
    acc = acc * 9.807#converts to m/s2
    amp = acc/(4*np.pi**2*peak_freq**2)

    return amp

def acc_spectrum(v):

    if type(v) is str:
        v1 = openfile(v)
    else:
        v1 = np.array(v.copy())

    spe, f = th.spectrum(v1, dt=1/freqwebdaq)

    return spe, f


def last_wdfile(ext=None):
    searchext = 'OBB*'
    if ext is not None:
        searchext = ext+'*'
    fl = sorted(glob.glob(basepathwebdaq+searchext))

    ffl = fl[-1]
    fname = ffl.split('/')[-1]
    print(fname)
    return fname


def last_N_wdfiles(N, ext=None):
    searchext = 'Piezo*'
    if ext is not None:
        searchext = ext+'*'
    fl = sorted(glob.glob(basepathwebdaq+searchext))
    file_list = []
    for j in range(N):
        ffl = fl[-j]
        fname = ffl.split('/')[-1]
        print(fname)
        file_list.append(fname)
    return file_list


def _openwdd(fname, manual_acq_data_len = None):
    hdr = {}
    with open(fname,"rb") as wdf:
        rawData = wdf.read(564)
        hdr['version'] = int.from_bytes(rawData[0:4], 'little')
        hdr['size'] = int.from_bytes(rawData[4:8], 'little')
        hdr['nchannels']= int.from_bytes(rawData[8:12], 'little')
        hdr['scan_rate']= int.from_bytes(rawData[12:20], 'little')
        hdr['start_time']= int.from_bytes(rawData[20:28], 'little')
        hdr['timezone']= rawData[28:44].decode("utf-8")
        json_hdr_size =  int.from_bytes(rawData[560:564], 'little')
        jsonRaw = wdf.read(json_hdr_size)
        hdr['json_hdr']=json.loads(jsonRaw)
        try:
            ndata = hdr['json_hdr']['jobDescriptor']['acquisition']['stopTrigger']['sampleCount']
        except KeyError:
            ndata = manual_acq_data_len # used for manual acquisition
        data_it = struct.iter_unpack('<d', wdf.read(ndata*hdr['nchannels']*8)) #4 because double precision 64 bit\n",
        tmp = np.asarray(list(data_it), dtype='double')
        data=tmp.reshape(int(tmp.size/4), 4)
        data = data.T
    return data


