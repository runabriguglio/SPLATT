import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits

import os
import glob
import subprocess


def read_buffer_data(TN:str = None):

    ip = _get_local_ip()
    if '193.206.155.43' in ip: # M4WS
        os.system('buffsync')
        SPLATT_BUFFER_FOLDER = '/mnt/jumbo/SPLATT/Buffer/'

    elif '193.206.155.220' in ip: # SPLATTWS
        SPLATT_BUFFER_FOLDER = '/home/labot/Desktop/Data/SPLATT/Buffer/'

    else:
        raise Exception(f"Local ip {ip} might not have access to buffer folders")

    freq = 1818 # [Hz]

    if TN is not None:
        where = os.path.join(SPLATT_BUFFER_FOLDER,TN)
    else:
        buffer_folder_list = sorted(glob.glob(SPLATT_BUFFER_FOLDER+'2025*'))
        where = buffer_folder_list[-1].split('/')[-1]

    dec = read_fits(where,'decimation.fits')
    if dec is None:
        raise FileNotFoundError('The TN does not seem to contain any decimation.fits file')

    dataR1 = read_fits(where,'dataR1.fits')
    dataR2 = read_fits(where,'dataR2.fits')
    dataW1 = read_fits(where,'dataW1.fits')
    dataW2 = read_fits(where,'dataW2.fits')
    
    startPosCmd = read_fits(where,'start_sabu16_position.fits')
    startCurCmd = read_fits(where,'start_sabi16_force.fits')

    print(np.shape(dataR1))
    data_len = np.shape(dataR1)[-1]
    dt = 1./freq*(dec+1.)
    time_vec = np.arange(data_len)*dt

    data = []
    data_addr = []
    data.append(dataR1)
    data_addr.append(_read_sab_address(where,'addrR1.fits'))
    if dataR2 is not None:
        data_addr.append(_read_sab_address(where,'addrR2.fits'))
        data.append(dataR2)
    if dataW1 is not None:
        data_addr.append(_read_sab_address(where,'addrW1.fits'))
        data.append(dataW1)
    if dataW2 is not None:
        data_addr.append(_read_sab_address(where,'addrW2.fits'))
        data.append(dataW2)
    if startPosCmd is not None:
        data_addr.append('start_sabu16_position')
        data.append(startPosCmd)
    if startCurCmd is not None:
        data_addr.append('start_sabi16_force')
        data.append(startCurCmd)

    return data, data_addr, time_vec


def spectral_analysis(signal, dt:float = 0):
    """ This is equivalent to:
        from M4/m4/mini_OTT import time_history as th
        th.spectrum(signal,dt)
        """

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

    # dt = (dec+1)/1818

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


def splatt_plot(values,min_val=None, max_val=None):
    coordAct = np.loadtxt('../SPLATT_Data/act_coords.txt')
    nActs = len(coordAct)

    # Perform matrix rotation to align with reference
    phi = 60./180*np.pi
    c=np.cos(phi)
    s=np.sin(phi)
    MatRot=[[c,-s],[s,c]]
    coordAct = MatRot@coordAct.T
    coordAct = coordAct.T

    # Set scatter plot variables
    Margin = 0.03
    markerSize = 800
    x = coordAct[:,0]
    y = coordAct[:,1]
    indices = np.arange(nActs)+1

    # Plot
    plt.figure()
    ax = plt.axes()
    ax.set_xlim(min(x)-Margin,max(x)+Margin)
    ax.set_ylim(min(y)-Margin,max(y)+Margin)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.scatter(x, y, c=values, vmin = min_val, vmax = max_val, s=markerSize, edgecolor='k')
    plt.colorbar()

    # Write 'G' reference and actuator indices
    for i in range(nActs):
        plt.text(x[i]*2/3, y[i]+Margin*2/3, str(indices[i]))
    plt.text(x[15],y[15]*1.3,'G')
    plt.show()


def mirror_mesh(values):
    npix = 128

    X = int(2*npix)
    Y = int(2*npix)
    circ_mask = np.fromfunction(lambda i,j: np.sqrt((i-npix)**2+(j-npix)**2) > npix, [X,Y])

    IFF = np.loadtxt('../SPLATT_Data/iffs.txt')

    flat_img = np.zeros(np.size(circ_mask))
    flat_img[~circ_mask.flatten()] = IFF @ values
    img = np.reshape(flat_img, np.shape(circ_mask))
    masked_img = np.ma.masked_array(img, circ_mask)

    plt.figure()
    plt.imshow(masked_img, origin='lower')


def _read_sab_address(folder_path, file_name):

    raw_addr = read_fits(folder_path, file_name)
    int_addr = (raw_addr[0]).astype(int)

    addr = chr(int_addr[0])
    for i in range(len(int_addr)-1):
        addr += chr(int_addr[i+1])

    return addr


def read_fits(file_path:str, file_name:str):

    which = os.path.join(file_path,file_name)
    try:
        hdu = pyfits.open(which)
        read_data =np.array(hdu[0].data)
    except FileNotFoundError:
        read_data = None

    return read_data


def _get_local_ip():

    res = subprocess.run(['hostname','-I'], stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, text=True)
    ip = res.stdout.strip().split(' ')

    return ip





