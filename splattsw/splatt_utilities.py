import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits

import os
import glob
import subprocess


def read_buffer_data(TN:str = None):

    ip = _get_local_ip()
    if '193.206.155.43' in ip: # M4WS
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
        try:
            print('Synchronizing buffer folder ...')
            subprocess.run("rsync -av --include='*/' --include='*.fits' --exclude='*' --prune-empty-dirs labot@splatt:/home/labot/Desktop/Data/SPLATT/Buffer/ /mnt/jumbo/SPLATT/Buffer/",shell=True)
            dec = read_fits(where,'decimation.fits')
        except:
            raise FileNotFoundError('The TN does not seem to contain any decimation.fits file')

    dataR1 = read_fits(where,'dataR1.fits')
    print(np.shape(dataR1))

    #startPosCmd = read_fits(where,'start_sabu16_position.fits')
    #startCurCmd = read_fits(where,'start_sabi16_force.fits')

    data_len = np.shape(dataR1)[-1]
    dt = 1./freq*(dec+1.)
    time_vec = np.arange(data_len)*dt

    # Build a data dictionary:
    addrR1 = _read_sab_address(where,'addrR1.fits')
    addrR2 = _read_sab_address(where,'addrR2.fits')
    addrW1 = _read_sab_address(where,'addrW1.fits')
    addrW2 = _read_sab_address(where,'addrW2.fits')
    data = {addrR1 : dataR1}
    if addrR2 is not None:
        dataR2 = read_fits(where,'dataR2.fits')
        data[addrR2] = dataR2
    if addrW1 is not None:
        dataW1 = read_fits(where,'dataW1.fits')
        data[addrW1] = dataW1
    if addrW2 is not None:
        dataW2 = read_fits(where,'dataW2.fits')
        data[addrW2] = dataW2

    return data, time_vec


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

    try:
        raw_addr = read_fits(folder_path, file_name)
        int_addr = (raw_addr[0]).astype(int)

        addr = chr(int_addr[0])
        for i in range(len(int_addr)-1):
            addr += chr(int_addr[i+1])
    except FileNotFoundError:
        addr = None

    return addr


def read_fits(file_path:str, file_name:str):
    
    try:
        which = os.path.join(file_path,file_name)
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





