from M4.m4.devices.i4d import I4D
from M4.m4.configuration.ott_parameters import Interferometer
interf = I4D(Interferometer.i4d_IP, Interferometer.i4d_port) 
from M4.m4.ground.timestamp import Timestamp
import os

basepath = 'D:/SPLATT/'
storefolder =  '/mnt/jumbo/SPLATT/OPTData/'
mountpoint4d = '/home/labot/4d/SPLATT/H5/'
freqfile = 'AcqInfo.txt'
#analysis of interferometer data
imgfile='/home/labot/testdir/999.4D'

def capture(nframes):
    
    tn = Timestamp.now()
    fold = basepath+'RAW/'+tn+'/raw/'
    interf.burstFramesToSpecificDirectory(fold, nframes)
    return tn

def produce(tn):

    rawFramesDirectory = basepath+'RAW/'+tn+'/raw/'
    dest = basepath+'H5/'+tn+'/hdf5/'
    interf.convertRawFramesInDirectoryToMeasurementsInDestinationDirectory(dest, rawFramesDirectory)

def frames_transfer(tn):
    print('Sync with 4D data folder')
    os.system('rsync -a '+mountpoint4d+tn+' '+storefolder)

def save_acqdata(tn,freq4d,piV, piFreq):
    logfile = storefolder+tn+'/'+freqfile
    f=open(logfile, 'w')
    dataw = ('%f\n%f\n%f\n\n#Freq 4d, PI command, PI Freq.' %(freq4d,piV,piFreq))
    f.write(dataw)
    f.close()

def read_acqdata(tn):
    print('Format: freq4d,PIcomm, PIfreq')
    logfile = storefolder+tn+'/'+freqfile
    fileexist = os.path.isfile(logfile)
    w = (0,0,0)
    if fileexist ==1:
        f=open(logfile,'r')
        w = []
        for i in range(3):
            s = f.readline()
            s = float(s)
            w.append(s)
        f.close()
    return w

def create_h5folder(tn):
    dest = storefolder+tn+'/'
    os.mkdir(dest)
