import numpy as np
from matplotlib import pyplot as plt

import os
import glob
import json
import struct
from astropy.io import fits as pyfits

from scipy.signal import lfilter

# WebDAQ variables
freqwebdaq = 1651.6129 #Hz; minimum sampling frequency
basepathwebdaq = '/mnt/jumbo/SPLATT/WebDaqData/'
#ftpwebdacq = '/home/ftpuser/ftp/' # old files
ftpwebdacq = '/home/ftpuser/ftp/files/' # new files, from 2024 Q4 onwards
webdext = 'OBB'

def wdsync():
    os.system('rsync -av -q -i '+ftpwebdacq+' '+basepathwebdaq)

def savefile(wdfname, tn):
    data = openfile(wdfname)
    pyfits.writeto(basepathwebdaq+tn+'.fits', data)


def openfile(filename, data_len = None):

    filext = filename.split('.')[-1]

    if filext == 'wdd':
        file_path = os.path.join(basepathwebdaq, filename)
        data = _openwdd(file_path, data_len)

    elif filext == 'fits':
        file_path = os.path.join(basepathwebdaq, filename)
        with pyfits.open(file_path) as hdul:
            data = hdul[0].data
    else:
        file_path = os.path.join(basepathwebdaq, filename + '.fits')
        with pyfits.open(file_path) as hdul:
            data = hdul[0].data

    return data

def plot_data(data, ch_ids = None, freq = freqwebdaq, title_str = None):
    data_size = np.shape(data)

    if ch_ids is None:
        ch_ids = np.arange(data_size[0])

    t_vec = np.arange(0,np.max(data_size))*1/freq

    plt.figure()
    if title_str is not None:
        plt.title(title_str)
    for k,ch in enumerate(ch_ids):

        if len(data_size) == 1:
            vec = data
        else:
            vec = data[int(ch)]

        plt.plot(t_vec,vec,label='Channel ' + str(ch))

    plt.legend()
    plt.xlabel('Time [s]')
    plt.ylabel('Acceleration [g]')
    plt.grid('on')
    plt.axis('tight')
    
def integrate_peak_freq(spe, f_vec, bound = None, bins = 2):

    _, peak_freq, peak_id = find_peak_freq(spe, f_vec, bound)

    N = len(peak_id)
    if N > 1:
        displ = np.zeros(N)
        for i in range(N):
            displ = acc_integrate(spe[i], peak_freq[i], peak_id[i], delta_peak = bins)
    else:
        displ = acc_integrate(spe, peak_freq, peak_id, delta_peak = bins)

    return displ, peak_freq


def find_peak_freq(spe, freq_vec, bound = None):

    freq_len = len(freq_vec)
    idf = np.arange(freq_len)

    if bound is not None:
        idf = np.where(np.logical_and(freq_vec>= bound[0], freq_vec<=bound[1]))

    size_spe = np.shape(spe)
    nChan = size_spe[0]

    if freq_len is not size_spe[1]:
        spe = spe.T
        nChan = size_spe[0]

    peak_id = np.zeros(nChan)
    peak_val = np.zeros(nChan)
    peak_freq = np.zeros(nChan)

    for j in range(nChan):
        v = spe[j]
        peak = np.max(v[idf])
        peak_val[j] = peak*np.sqrt(2) # re-normalize amplitude
        peak_id = np.argmax(v[idf])
        idf1 = idf[0]
        peak_freq[j] = freq_vec[idf1[peak_id]]

    return peak_val, peak_freq, peak_id


def acc_integrate(spe, peak_freq, peak_id, delta_peak = 3):
    #acc is the acceleration value in g
    min_id = np.max((int(0),peak_id-delta_peak))
    max_id = np.min(len(spe),peak_id+delta_peak+1)
    spe_peak = spe[min_id:max_id]
    acc = np.sum(spe_peak)
    acc = acc * 9.807#converts to m/s2
    amp = acc/(4*np.pi**2*peak_freq**2)

    return amp

def acc_spectrum(v):

    if type(v) is str:
        v1 = openfile(v)
    else:
        v1 = np.array(v.copy())

    spe, f = get_spectrum(v1, dt=1/freqwebdaq)

    return spe, f


def last_wdfile(N:int=1, ext=None):
    searchext = webdext+'*' 
    if ext is not None:
        searchext = ext+'*'
    fl = sorted(glob.glob(basepathwebdaq+searchext))

    file_list = []
    for j in range(N):
        ffl = fl[-j-1]
        fname = ffl.split('/')[-1]
        print(fname)
        file_list.append(fname)

    if N == 1:
        file_list = file_list[0]


    return file_list


def get_spectrum(signal, dt=1, phase = False):
    spe = np.fft.rfft(signal, norm="ortho", axis=-1)
    nn = np.sqrt(spe.shape[-1])
    spe_mod = (np.abs(spe)) / nn
    freq = np.fft.rfftfreq(signal.shape[-1], d=dt)

    if len(np.shape(spe))>1:
        spe_mod[:,0] = 0
    else:
        spe_mod[0] = 0

    if phase is True:
        spe_phi = np.angle(spe)
        return spe_mod, freq, spe_phi
    else:
        return spe_mod, freq


def damp_sinusoid_coeffs(signal):

    spe = np.fft.rfft(signal, norm="ortho", axis=-1)
    spe = (np.abs(spe)) / np.sqrt(spe.shape[-1])
    spe[0] = 0  

    k0 = np.argmax(np.abs(spe))
    N = len(signal)
    n = np.arange(N)
    expk = lambda k:  np.exp(-1j*2*np.pi/N*(k)*n)
    Xp = np.dot(signal, expk(k0+0.5))
    Xm = np.dot(signal, expk(k0-0.5))

    sumX = Xp+Xm
    diffX = Xp-Xm

    est_phi = 0.5*np.real(sumX/diffX)
    est_damp = np.pi/N*np.imag(sumX/diffX)
    est_freq = 2*np.pi/N*(k0+est_phi)

    c1 = 2*np.cos(est_freq)
    c2 = -np.exp(-2*est_damp)

    # Filter step
    b = np.array([1])
    a = np.array([1, -c1, -c2])

    V = lfilter(b,a,signal)

    delta = np.zeros(N)
    delta[0] = 1
    U = lfilter(b,a,delta)

    # prediction matrix
    PM = np.zeros([N-2,4])
    PM[:,0] = V[1:N-1]
    PM[:,1] = V[0:N-2]
    PM[:,2] = U[2:N]
    PM[:,3] = U[1:N-1]

    b = V[2:]

    c = np.linalg.pinv(PM) @ b

    omega = np.arccos(c[0]/(2*np.sqrt(-c[1])))
    eta = -np.log(-c[1])/2

    exp1 = (1-np.exp(-eta))/(1-np.exp(-eta*N))
    exp2 = (1-np.exp(-(1j*2*omega+eta)*N))/(1-np.exp(-(1j*2*omega+eta)))
    sigsum = np.dot(signal,np.exp(-1j*omega*n))
    AA = lambda Ast:  exp1 * (sigsum-Ast*exp2)

    A_est = 0
    tol = 1e-12
    err = tol + 1
    max_it = 10
    it = 0
    while np.logical_and(it < max_it, err > tol):
        A_old = A_est
        A_est = AA(A_est)
        err = np.abs(A_old-A_est)
        it+=1

    if err > tol:
        raise ValueError(f'The iteration has not treach convergence after {max_it} iterations: the residual {err} is still greater than the tolerance {tol}')

    amp = 2*np.sqrt(np.imag(A_est)**2+np.real(A_est)**2)
    phi = np.arctan(np.imag(A_est)/np.real(A_est))

    return amp, eta, omega, phi


def find_peaks(signal, Npeaks:int=2, min_sep:int=10):

    abs_sig = np.abs(signal)
    sorted_ids = np.argsort(abs_sig)

    ctr = -1
    peak_ids = np.zeros(Npeaks,dtype=int)
    peak_ids[0] = sorted_ids[ctr]

    for ii in range(Npeaks-1):
        ctr -= 1
        peak_id = sorted_ids[ctr]

        while np.min(np.abs(peak_id-peak_ids)) < min_sep:
            ctr -= 1
            peak_id = sorted_ids[ctr]

        peak_ids[ii+1] = peak_id

    return peak_ids



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


