#from SPLATT.splattsw import splatt_analysis as sp
from M4.m4.ground.timestamp import Timestamp
import os
import numpy as np
basepath = '/mnt/jumbo/SPLATT/'
logfile = basepath+(Timestamp.now())[0:8]+'.log'
sep = '\t'

def write_header():
    fexist= os.path.isfile(logfile)

    if fexist == False:
        f=open(logfile,'a')
        w = ''
        s = ['Tracknum 4d', 'WebDAQ file', 'LATT Buffer', 'Freq4d','PI Comm', 'PI Freq','Peak OBB[um]','Peak Stand [um]','Peak TT [um]','PhaseIn', 'PhaseOBB', 'PhaseStand', 'PhaseTT']
        for i in s:
            w = ('%s%s%s' %(w,i,sep))
        w = w+'\n'
        f.write(w)
        f.close()

def log_data(datafile, dataval):
    write_header()
    f=open(logfile,'a')
    sep = '\t'
    nn = np.size(datafile)
    if (nn != 3):
        datafile.append('None')

    s = ''
    for i in datafile:
        s = ('%s%s%s' %(s,i,sep))
    
    #s = ('%s\t%s\t%s\t' %(datafile[0],datafile[1],datafile[2]))
    for i in dataval:
        s = ('%s%f%s' %(s,i,sep))
    s = s+'\n'
    f.write(s)
    f.close()
    print('Info logged to '+logfile)



