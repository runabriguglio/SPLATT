from SPLATT.splattsw import splatt_analysis as sp
from SPLATT.splattsw.devices.webDAQ import WebDAQ as wbdq
from SPLATT.splattsw.devices import wavegenerators as wg

import matplotlib.pyplot as plt
import numpy as np
import time


# Connect to WebDAQ
webdaq = wbdq()
#webdaq.connect()
freq_vec = np.arange(10,40+3,3)

time_for_acquisition = 9.
file_string_initials = 'OBB'

def find_peak(v, freq=None, bound=None):
    #requires a monodimensional array
    idf = range(len(v))
    if freq is not None and bound is not None:
        idf =     np.where(np.logical_and(freq>= bound[0], freq<=bound[1]))
    peak = max(v[idf])
    peakid = np.argmax(v[idf])
    peakfreq = 0
    if freq is not None:
        idf1 = idf[0]
        peakfreq = freq[idf1[peakid]]

    return peak, peakfreq, idf1[peakid]

def acc_integrate(spe, peak_freq, peak_id, delta_peak = 3):
    #acc is the acceleration value in g
    #print('Is the acceleration provided as [g]?')
    spe_peak = spe[(peak_id-delta_peak):(peak_id+delta_peak+1)]
    acc = np.sum(spe_peak)
    acc = acc * 9.807#converts to m/s2
    amp = acc/(4*np.pi**2*peak_freq**2)

    return amp


def analyse_wdfile(wdfile, exc_freq, doplot = False, ch = 0):
    #accelerometer analysis
    #spectrum
    spe_4ch, f = sp.acc_spectrum(wdfile)
    spe = spe_4ch[ch,:]*np.sqrt(2) # re-multiply by sqrt(2) to have amp=oscillation_amp

    if doplot == True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(f,spe)
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel("Amplitude [g]")
        fig.show()

    peak, peak_f, peakId = find_peak(spe,freq=f,bound=[exc_freq-2.,exc_freq+2.])#sp.find_peak(spe,freq=f,bound=[1.,f[-1]])
    peak_d = acc_integrate(spe,peak_f,peakId) #sp.acc_integrate(peak,peak_f)

    return peak_d, peak_f


def test_single_freq(freq_val,amp=1):
    # Connect to WebDAQ
    # webdaq.connect()

    wg.set_wave1(amp,0,freq_val,'SIN')
    time.sleep(2) # wait for piezo command
    webdaq.start_schedule()
    time.sleep(time_for_acquisition) # wait for acquisition to end

    wg.clear_wave1()

    #read file
    sp.wdsync() # does os.system('rsync -av '+ftpwebdacq+' '+basepathwebdaq)
    wdfile = sp.lastwdfile(file_string_initials)

    return wdfile


def start_cycle_on_freq(freqV = freq_vec):
    ampPI = 1
    scale_fact = ampPI
    file_list = []

    # Connect to WebDAQ
    webdaq.connect()

    # Select device ('Rigol_WaveGen' or 'RedPitaya')
    wg.update_device('Rigol_WaveGen')

    N = len(freqV)
    Nch = 4
    maxD = np.zeros([N,Nch])
    maxF = np.zeros([N,Nch])

    for i in range(N):
        freqPI = freqV[i]
        wdf=test_single_freq(freqPI,ampPI) #ampPi/ampGain
        for j in range(Nch):
            maxD[i,j], maxF[i,j] = analyse_wdfile(wdf,freqPI,ch=j)
        file_list.append(wdf)

    for j in range(Nch):
        maxD[:,j] = maxD[:,j]/scale_fact
        plt.figure()
        plt.plot(maxF[:,j], maxD[:,j])
        plt.scatter(maxF[:,j], maxD[:,j], marker='o', c='red', s=15)
        plt.xlabel("Peak frequency [Hz]")
        plt.ylabel("Peak oscillation [m]")
        plt.title("Channel "+str(j))
        plt.show()

    return file_list, maxD, maxF


def plot_time_data(wdfile,freq_WebDAQ = 1651.,Nch=4):
    data = sp.openfile(wdfile)

    data_mean = np.mean(data,axis=1)
    data_mean = np.repeat(data_mean,len(data[0]))
    data_mean = data_mean.reshape(np.shape(data))

    data_norm = data - data_mean

    dt = 1./freq_WebDAQ
    t_vec = np.arange(0,dt*len(data[0,:]),dt)

    plt.figure()
    for i in range(Nch):
        plt.plot(t_vec,data_norm[i])
    plt.xlabel('Time [s]')
    plt.ylabel('Acceleration [g]')
    plt.legend(['OBB','Stand','Piezo','Spring'])
    plt.show()

    return data_norm

def analyse_file_list(file_list,freq_vec,Nch=4):

    L = len(file_list)
    maxD = np.zeros([Nch,L])
    maxF = np.zeros([Nch,L])

    for k in range(len(file_list)):
        wdf = file_list[k]
        freq = freq_vec[k]
        for chan in range(Nch):
            maxD[chan,k],maxF[chan,k] = analyse_wdfile(wdf,freq,ch=chan)

    return maxD, maxF

# def recover_planar_oscillation(wdfile, freq):
#     '''Only works when the accelerometers are placed at 90°
#     from each other on the OBB'''
#
#     data = sp.openfile(wdfile)
#
#     data_mean = np.mean(data,axis=1)
#     data_mean = np.repeat(data_mean,len(data[0]))
#     data_mean = data_mean.reshape(np.shape(data))
#
#     dat = data - data_mean
#
#     T = 1./freq
#
#     freq_WebDAQ = 1651.
#     dt = 1./freq_WebDAQ
#
#     time_vec = np.arange(0,dt*len(data[0,:]),dt)
#
#     shift = int(np.round(T/4./dt)) # corresponds to a phase shift of T/4
#     print(shift)
#
#     t_vec = time_vec[:-shift]
#
#     sdat = np.zeros([4,len(data[0,:])-shift])
#     sdat[0,:] = dat[0,shift:]
#     sdat[1,:] = dat[1,:-shift]
#     sdat[2,:] = dat[2,shift:]
#     sdat[3,:] = dat[3,:-shift]
#
#     ddalpha = (sdat[0]-sdat[2])*9.807/0.20/2 # oscillation about axis
#     ddbeta = (sdat[1]-sdat[3])*9.807/0.20/2 # oscillation along axis
#     ddtheta = 9.807*((sdat[0]-sdat[2])/0.21/2 + (sdat[1]+sdat[3])/0.07/2)
#
#     plt.figure()
#     plt.plot(t_vec,ddalpha)
#     plt.plot(t_vec,ddbeta)
#     plt.plot(t_vec,ddtheta)
#     plt.legend(['x','y','axis'])
#     plt.xlabel('Time [s]')
#     plt.ylabel('Angular acceleration [s^-2]')
#     plt.show()
#
#     return ddtheta


# vec = np.zeros(10)
# vecF=np.zeros(10)
#
# for i in range(10):
#     data = sp.openfile(fl[i])
#     vec[i] = np.std(data[0])
#     spe, f = sp.acc_spectrum(fl[i])
#     peak, peak_f, pid = find_peak(spe[0],freq=f,bound=[1.,f[-1]])
#     vec[i]=vec[i]/(4*np.pi**2*peak_f**2)
#     vecF[i]=peak_f
#     print(peak_f)
#
#
#
# H = maxD[:,0]/maxD[:,1]
#
# plt.figure()
# plt.plot(maxF[:,0], 20.*np.log10(H))
# plt.xlabel("Frequency [Hz]")
# plt.ylabel("OBB/stand vibration [dB]")
# plt.title('Stand to OBB transfer function')
# plt.grid('on')
# plt.show()
#
# H = maxD[:,0]/maxD[:,2]
#
# plt.figure()
# plt.plot(maxF[:,0], 20.*np.log10(H))
# plt.xlabel("Frequency [Hz]")
# plt.ylabel("OBB/piezo vibration [dB]")
# plt.title('Piezo to OBB transfer function')
# plt.grid('on')
# plt.show()
#
# H = maxD[:,3]/maxD[:,2]
#
# plt.figure()
# plt.plot(maxF[:,0], 20.*np.log10(H))
# plt.xlabel("Frequency [Hz]")
# plt.ylabel("arm/piezo vibration [dB]")
# plt.title('Piezo to stand arm bottom transfer function')
# plt.grid('on')
# plt.show()
#
# plt.figure()
# for i in range(4):
#     plt.plot(maxF[:,i], maxD[:,i])
# plt.xlabel("Frequency [Hz]")
# plt.ylabel("Oscillation [m]")
# plt.grid('on')
# plt.show()
#





