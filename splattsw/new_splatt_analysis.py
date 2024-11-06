import numpy as np
from matplotlib import pyplot as plt
import os
import json
import struct

freqwebdaq = 1651.6129 #Hz; minimum sampling frequency
basepathwebdaq= '/mnt/jumbo/SPLATT/WebDaqData/'

def openfile(name):
    file_path = os.path.join(basepathwebdaq, name)
    data = _openwdd(file_path)
    return data

def plot_data(data, freq = freqwebdaq):

    data_size = np.shape(data)

    N_ch = data_size[0]
    N_tvec = data_size[1]

    t_vec = np.arange(0,N_tvec)*1/freq

    for i in range(N_ch):
        plt.figure()
        plt.plot(t_vec,data[i])
        plt.title('Channel '+str(i))
        plt.xlabel('Time [s]')
        plt.ylabel('Acceleration [g]')
        plt.show()


def _openwdd(fname):
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
        ndata = hdr['json_hdr']['jobDescriptor']['acquisition']['stopTrigger']['sampleCount']
        data_it = struct.iter_unpack('<d', wdf.read(ndata*hdr['nchannels']*8)) #4 because double precision 64 bit\n",
        tmp = np.asarray(list(data_it), dtype='double')
        data=tmp.reshape(int(tmp.size/4), 4)
        data = data.T
    return data


