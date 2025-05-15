import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits

from acceleration_analysis import get_spectrum

import os
import glob
import subprocess

def read_fits(file_path:str, file_name:str):
    
    try:
        which = os.path.join(file_path,file_name)
        with pyfits.open(which) as hdu:
            hdu.verify('fix')
            read_data =np.array(hdu[0].data)
    except FileNotFoundError:
        read_data = None

    return read_data


def buffsync():
    print('Synchronizing buffer folder ...')
    subprocess.run("rsync -av --include='*/' --include='*.fits' --exclude='*' --prune-empty-dirs labot@splatt:/home/labot/Desktop/Data/SPLATT/Buffer/ /mnt/jumbo/SPLATT/Buffer/", shell=True)


def read_buffer(TN:str = None):

    ip = _get_local_ip()
    if '193.206.155.43' in ip: # M4WS
        SPLATT_BUFFER_FOLDER = '/mnt/jumbo/SPLATT/Buffer/'

    elif '193.206.155.220' in ip: # SPLATTWS
        SPLATT_BUFFER_FOLDER = '/home/labot/Desktop/Data/SPLATT/Buffer/'

    else:
        raise Exception(f"Local ip {ip} might not have access to buffer folders")

    freq = 1/550e-6 # [Hz]

    if TN is not None:
        where = os.path.join(SPLATT_BUFFER_FOLDER,TN)
    else:
        buffer_folder_list = sorted(glob.glob(SPLATT_BUFFER_FOLDER+'2025*'))
        where = buffer_folder_list[-1].split('/')[-1]

    dec = read_fits(where,'decimation.fits')
    if dec is None:
        raise FileNotFoundError('The TN does not seem to contain any decimation.fits file')

    dataR1 = read_fits(where,'dataR1.fits')
    print(np.shape(dataR1))

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


def analyse_buffer(TN:str = None):

    data_dict, tvec = read_buffer(TN)
    ZM = _get_splatt_zernike_matrix()
    tvec = tvec.T
    dt = tvec[1]-tvec[0]

    data = data_dict.copy()

    for key in data_dict.keys():
        val = data_dict[key].copy()
        if np.shape(val)[1] < np.shape(val)[0]:
            val = val.T
        spe, fvec = get_spectrum(val, dt)
        data[key + '_spectrum'] = spe
        zdata = ZM.T @ val
        dzdata = zdata - np.reshape(zdata[:,0],[np.shape(zdata)[0],1])
        data[key + '_zernike'] = dzdata
        spe_z, _ = get_spectrum(dzdata, dt)
        data[key + '_zernike_spectrum'] = spe_z

    # Convert everything in a format friendlier to matplotlib.pyplot
    for key in data.keys():
        data[key] = data[key].T

    tvec = np.repeat(tvec,19,axis=1)

    return data, tvec, fvec


def splatt_plot(values, min_val=None, max_val=None):
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
    except: #FileNotFoundError:
        addr = None

    return addr


def _get_local_ip():

    res = subprocess.run(['hostname','-I'], stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, text=True)
    ip = res.stdout.strip().split(' ')

    return ip


def _get_splatt_zernike_matrix():

    coordAct = np.loadtxt('../SPLATT_Data/act_coords.txt')
    nActs = len(coordAct)

    ZM = np.zeros([nActs, nActs])

    for i in range(nActs):
        ZM[:,i] = _project_zernike(i+1, coordAct)

    return ZM


def _project_zernike(noll_number:int, coords):
    """
    Project the Zernike polynomials identified by the Noll number in input
    on the given coordinates.
    The polynomials are computed on the largest circle in the coordinates

    Parameters
    ----------
    noll_number : int
        Noll index of the desired Zernike polynomial.
    coords : [Npoints,2]
        Coordinates of the desired projection.

    Returns
    -------
    zern_data : [Npoints]
        Flattened array of the Zernike shape on the coordinates

    """
    if noll_number < 1:
        raise ValueError("Noll index must be equal to or greater than 1")

    # Image dimensions
    X = coords[:,0]
    Y = coords[:,1]

    # Determine circle radius on to which define the Zernike
    r = np.max(np.sqrt(X**2+Y**2))/2

    # Conversion to polar coordinates on circle of radius r 
    phi = lambda i,j: np.arctan2((j-Y/2.)/r,(i-X/2.)/r)
    rho = lambda i,j: np.sqrt(((j-Y/2.)/r)**2+((i-X/2.)/r)**2)
            
    zmode = lambda i,j: _zernikel(noll_number, rho(i,j), phi(i,j))
    mode = zmode(X,Y)

    # Normalization of the masked data: null mean and unit STD
    if noll_number > 1:
        norm_mode = (mode - np.mean(mode))/np.std(mode)
    else:
        norm_mode = mode

    return norm_mode
    

########### from M4SW ############

def _zernikel(j, rho, phi):
    """
    Calculate Zernike polynomial with Noll coordinate j given a grid of radial
    coordinates rho and azimuthal coordinates phi.

    >>> zernikel(0, 0.12345, 0.231)
    1.0
    >>> zernikel(1, 0.12345, 0.231)
    0.028264010304937772
    >>> zernikel(6, 0.12345, 0.231)
    0.0012019069816780774
    """
    m, n = _j2mn_noll(j)
    
    return _zernike(m, n, rho, phi)


def _zernike(m, n, rho, phi):
    """
    Calculate Zernike polynomial (m, n) given a grid of radial
    coordinates rho and azimuthal coordinates phi.

    >>> zernike(3,5, 0.12345, 1.0)
    0.0073082282475042991
    >>> zernike(1, 3, 0.333, 5.0)
    -0.15749545445076085
    """
    if (m > 0): return _zernike_rad(m, n, rho) * np.cos(m * phi)
    if (m < 0): return _zernike_rad(-m, n, rho) * np.sin(-m * phi)
    return _zernike_rad(0, n, rho)


from math import factorial as fac
def _zernike_rad(m, n, rho):
    """
    Calculate the radial component of Zernike polynomial (m, n)
    given a grid of radial coordinates rho.

    >>> zernike_rad(3, 3, 0.333)
    0.036926037000000009
    >>> zernike_rad(1, 3, 0.333)
    -0.55522188900000002
    >>> zernike_rad(3, 5, 0.12345)
    -0.007382104685237683
    """

    if (n < 0 or m < 0 or abs(m) > n):
        raise ValueError

    if ((n-m) % 2):
        return rho*0.0

    pre_fac = lambda k: (-1.0)**k * fac(n-k) / ( fac(k) * fac( int((n+m)/2.0 - k) ) * fac( int((n-m)/2.0 - k) ) )

    return sum(pre_fac(k) * rho**(n-2.0*k) for k in range((n-m)//2+1))


def _j2mn_noll(j):
    """
    Find the [n,m] list giving the radial order n and azimuthal order
    of the Zernike polynomial of Noll index j.

    Parameters:
        j (int): The Noll index for Zernike polynomials

    Returns:
        list: n, m values
    """
    n = int((-1.+np.sqrt(8*(j-1)+1))/2.)
    p = (j-(n*(n+1))/2.)
    k = n%2
    m = int((p+k)/2.)*2 - k

    if m!=0:
        if j%2==0:
            s=1
        else:
            s=-1
        m *= s

    return [m, n]





