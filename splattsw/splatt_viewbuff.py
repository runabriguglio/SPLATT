from SPLATT.splattsw import splatt_analysis as sp
from astropy.io import fits as pyfits
from M4.m4.mOTT_analysis import timehistory as th
from importlib  import reload #for reload
import scipy.io

mat = scipy.io.loadmat('file.mat')
