import numpy as np
import sys

sys.path.insert(0,'../..')

from M4.m4.ground import timestamp as tt
from SPLATT.splattsw.devices import moxa_io as mx

def log_temperature():
    """ Function to log the temperature"""

    tn=tt.Timestamp()
    time = tn.now()
    pt1 = mx.moxa_ai('PT1')
    temp = pt1.read()

    info_line = f'{time}     {str(temp[0:3])}\n'

    fname='/home/labot/git/SPLATT/splattsw/temp_reads.txt'
    fff=open(fname, 'a+')
    fff.write(info_line)
    fff.close()
