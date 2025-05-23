import os
import glob
import numpy as np
#import jdcal
from astropy.io import fits as pyfits
#from m4.ground import read_data
#from m4.ground import zernike
#from m4.ground.read_data import InterferometerConverter
#ic = InterferometerConverter()

from matplotlib import pyplot as plt

a= '/mnt/data/M4/Data/M4Data/OPTData/'



class TimeHist():
    '''
    Class to 
    HOW TO USE IT::

        
    '''

    def __init__(self, tn):
        """The constructor """
        self.tracknum = tn
         
        self._fold  = findTracknum(tn)
        self._path = a+ self._fold
        self._list = fileList(tn)
        
    def frame(self, id):
        
        return frame(id, self._list)
        
    def averageFrames(start, stop):
    
        return averageFrames(start, stop, self._list)

    

    #@staticmethod
    #def _storageFolder():
    #    """ Creates the path where to save data"""
    #    return fold_name.OPDSERIES



def findTracknum(tn):
    '''
    Parameters
    ----------
    tn: string
        tracking number to be searched in the data folder

    Returns
    -------
    result: the specific data folder where the tracknum is found
    '''

    #a= '/mnt/data/M4/Data/M4Data/OPTData/'
    lsdir = os.listdir(a)
    for i in lsdir:
        b = a+i
        z = os.listdir(b)
        check = False
        for j in z:
            check = (j == tn)
            if check == True:
                result = i
                return result

def fileList(tn, fold=None):
    '''
    Parameters
    ----------
    tn: string
        tracking number where to search for the images file list

    Returns
    -------
    lsdir: the list of image files
    fold1: the path where the file is found
    '''
    if fold is not None:
        name = '*.4D'
        addfold ='/hdf5/'
    else:
        
        fold = findTracknum(tn)
        addfold = '/'
        name = '20*'
        if fold == 'OPDImages':
            addfold = '/hdf5/'
            name = 'img*'
            
    fold1 = fold+'/'+tn+addfold   #to be re-checked at OTT!! 
    lsdir = sorted(glob.glob(fold1+name), key=lambda x: int(os.path.basename(x).split('/')[-1].split('.')[0]))
    #lsdir = lsdir[0]
     
    return lsdir

def read_phasemap(filename, thefold = None):

    #if thefold is not None:
    #    the = filename.split(thefold)[1]

    thetype = filename.split('.')[1]
    if thetype == 'fits':
        hduList = pyfits.open(filename)
        img = hduList[0].data
        mask = np.zeros(img.shape, dtype=np.bool)
        mask[np.where(img == img.max())] = True
        img = np.ma.masked_array(img, mask=mask)
       
        
        hduList.close()
        
    if thetype == '4D':
        print('4D')
        img = ic.fromNew4D(filename)

    if thetype == 'h5':
        img = ic.from4D(filename)
    
    return img

def averageFrames(first, last, fileList, thresh=None):
    '''
    Parameters
    ----------
    first: first item 
        tracking number where to search for the images file list

    Returns
    -------
    lsdir: the list of image files
    fold1: the path where the file is found
    '''
    imcube = cubeFromList(fileList[first:last+1])
    if thresh is None:
        aveimg = np.ma.mean(imcube, axis=0)
        
    else:
        img = imcube[0].data*0
        mmask = imcube[0].mask
        mysize = imcube[0].compressed()
        nn = 0
        for i in imcube:
            if i.data.compressed.size > 1:
                nn +=1
                img += i.data
                mmask = np.ma.mask_or(i.mask, mmask)
            
        img = img/nn
        image = np.ma.masked_array(img, mask=mmask)
        
    return aveimg

def removeZernike(ima, modes=np.array([1,2,3,4])):

        coeff, mat = zernike.zernikeFit(ima, modes)
        surf = zernike.zernikeSurface(ima, coeff, mat)
        new_ima = ima-surf
        return new_ima
        
def zernikePlot(mylist, modes=np.array(range(1,11))):
    mytype = type(mylist)
    if mytype is list:
        imgcube = cubeFromList(mylist)
        
    if mytype is np.ma.core.MaskedArray:
        imgcube = mylist
        
    zlist = []
    for i in range(len(imgcube)):
        coeff, mat = zernike.zernikeFit(imgcube[i], modes)
        zlist.append(coeff)
        
    zcoeff = np.array(zlist)
    zcoeff = zcoeff.T
    return zcoeff
    
def runningDiff(tn, gap=2):
    llist =fileList(tn)
    nfile = len(llist)
    npoints = nfile/gap-2
    slist=[]
    for i in range(0,npoints):
        #print(i*gap)
        #print(i*gap+1)
        q0 = frame(i*gap,llist)
        q1 = frame(i*gap+1,llist)
        diff = q1-q0
        diff = removeZernike(diff)
        slist.append(diff.std())
    svec = np.array(slist)
    return svec
    
def frame(id, mylist):

    mytype = type(mylist)
    if mytype is list:
        img = read_phasemap(mylist[id])
        
    if mytype is np.ma.core.MaskedArray:
        img = mylist[id]
    
    return img
    
def spectrum(signal, dt=1, show=None):
    
    # Spectrum
    spe  = np.fft.rfft(signal, norm='ortho') 

    # Normalization
    nn = 1
    nsig = signal.shape
    if np.size(nsig) > 1:
        nn   = spe.shape[1]
        
    spe  = (np.abs(spe)) / np.sqrt(nn)

    # Remove first element (zero frequency)
    if np.size(nsig) ==1:
        spe[0] = 0
    else:
        spe[:,0] = 0
    
    # Frequency vector
    freq = np.fft.rfftfreq(nsig[0], d=dt)
        
    if show is not None:
        for i in range(0,nn):
            plt.plot(freq, spe[i,:])
            
    return spe, freq
    
        
def cubeFromList(fileList):
    image_list = []
    for i in fileList:
        ima = read_phasemap(i)
        image_list.append(ima)

    image_list = np.ma.masked_array(image_list)
    return image_list


def timevec(tn):
    fold = findTracknum(tn)
    flist = fileList(tn)
    nfile = len(flist)
    if fold == 'OPDImages':
        tspace = 1./28.57
        timevec = range(nfile)*tspace
    
    
    if fold == 'OPD_series':
        timevec = []
        for i in flist:
            pp = i.split('.')[0]
            tni = pp.split('/')[-1]
            y=tni[0:4]
            mo = tni[4:6]
            d = tni[6:8]
            h = float(tni[9:11])
            mi = float(tni[11:13])
            s = float(tni[13:15])
            jdi=sum(jdcal.gcal2jd(y, mo, d))+h/24+mi/1440+s/86400
            timevec.append(jdi)
        timevec=np.array(timevec)
    
        
    return timevec
        
def track2jd(tni):
    y, mo, d, h, mi, s = track2date(tni)
    jdi=sum(jdcal.gcal2jd(y, mo, d))+h/24+mi/1440+s/86400  
    return jdi
    
def track2date(tni):
    y=tni[0:4]
    mo = tni[4:6]
    d = tni[6:8]
    h = float(tni[9:11])
    mi = float(tni[11:13])
    s = float(tni[13:15])
    return y, mo, d, h, mi, s
    
def runningMean(vec, npoints):
    
    return np.convolve(vec, np.ones(npoints), 'valid') / npoints       
        
        
        


