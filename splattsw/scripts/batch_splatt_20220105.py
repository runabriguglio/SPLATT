import os
import glob
import numpy as np
#import jdcal
from astropy.io import fits as pyfits
from SPLATT.splattsw.devices import webDAQ as wd
from M4.m4.mOTT_analysis import timehistory as th
from matplotlib.pyplot import *
from SPLATT.splattsw import splatt_log as slog
from SPLATT.splattsw import splatt4dmeas as comm4d
from SPLATT.splattsw import splatt_analysis as sp

#measurement of astigmatism to test if shell is detaching at 100 Hz
tn0='20211229_152836' #shell RIP, 40 Hz
tn1 = '20211229_154247' # shell rip, 100 Hz

tt0=sp.tiltvec(tn0)
tt1=sp.
