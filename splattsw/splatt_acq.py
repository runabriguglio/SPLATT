
from splattsw import acceleration_analysis as sp
from splattsw.devices.wavegenerators import WaveGenerator as wg
from splattsw.devices.webDAQ import WebDAQ as wdq

from opticalib.devices.interferometer import PhaseCam
from splattsw.devices.moxa_io import Moxa_ai0
from splattsw.devices.deformable_mirror import SPLATTEngine

from threading import Thread
import time
import numpy as np
from datetime import datetime
from astropy.io import fits as pyfits
import os

class Acquisition():

    def __init__(self):

        self.config_path = '/mnt/libero/SPLATTData/TestConfig'

        # Define interferometer
        self.interf = None
        try:
            self.interf=PhaseCam('4020')
        except:
            print('Interferometer not found')

        # Accelerometers
        try:
            self.webdaq=wdq()
            self.webdaq.connect()
        except:
            print('WebDAQ not found')

        # Signal generator
        self.wavegen = wg()

        # SPLATT dm
        self.dm = None
        try:
            self.dm = SPLATTEngine()
            self.eng = self.dm._eng
        except:
            print('SPLATT dm not found')

        # MOXA
        self.mx = None
        try:
            self.mx = Moxa_ai0()
        except:
            print('Moxa device not found')



    def acq_sweep(self, fmin = 30, fmax = 110, duration = 11, ampPI = 2, 
                  nframes:int = 2250, chPI:int = 1, produce:bool = False):

        # Generate new tn
        tn = self._generate_tn()

        # Prepare sweep parameters
        self.wavegen.set_wave(ch=chPI,ampl=ampPI,offs=0,freq=fmin,wave_form='SIN')
        time.sleep(4) # wait for steady state
        self.wavegen.set_sweep(chPI,fmin,fmax,duration,amp=ampPI)

        # Start acquisition and sweep
        self.webdaq.start_schedule()
        self.wavegen.trigger_sweep(chPI)
        time.sleep(0.5)
        self.interf.capture(nframes,tn)

        self.webdaq.stop_schedule_when_job_ends(job_id = 0) # Wait for webDAQ acquisition
        self.sync_and_save_last_wdfile(tn)

        if self.mx is not None:
            self.mx.save_pressure(self.config_path,tn)

        if self.dm is not None:
            self.dm.save_state(self.config_path,tn)

        if produce:
            self.interf.produce(tn)

        return tn


    def acq_freq(self, freqPI, freq4D, ampPI=1, nframes:int = 1000, produce:bool = False, buffer_dec:int = None):

        if ampPI > 1:
            raise ValueError(f'Amplitude {ampPI} is greater than 1, this might be outside interferometer capture range when exciting around the 20 Hz frequency')

        # Generate new tn
        tn = self._generate_tn()

        # Start acquisition and pulse trigger
        self.wavegen.set_wave(ch=1, ampl=ampPI, offs=0, freq=freqPI, wave_form='SIN')
        self.wavegen.trigg4D(freq4D)

        # Start buffer acquisition
        if buffer_dec is not None:
            self.eng.send(f'clear opts; opts.dec = {buffer_dec}; opts.sampleNR = 256; opts.save2mat = 0; opts.save2fits = 1; opts.tn = {tn}')
            def start_buffer():
                self.eng.send("splattAcqBufInt({'sabi32_Distance','sabi32_pidCoilOut'},opts)")
            buffer_thread = Thread(target = start_buffer)
            buffer_thread.start()
        
        # Wait for transient before starting WebDAQ
        time.sleep(4) # wait for transient
        self.webdaq.start_schedule()

        # Capture frames
        if nframes > 0:
            self.interf.setTriggerMode(True)
            self.interf.capture(nframes, tn)
        
        # Wait for webDAQ acquisition to end
        self.webdaq.stop_schedule_when_job_ends(job_id = 0)

        # Acquisition end
        self.wavegen.wave_off(1)
        self.wavegen.wave_off(2)
        
        # Post-processing
        self.sync_and_save_last_wdfile(tn)

        if self.mx is not None:
            self.mx.save_pressure(self.config_path,tn)
        if nframes > 0:
            self.interf.setTriggerMode(False)
            if produce:
                self.interf.produce(tn)
            dirpath = os.path.join(self.config_path,tn)
            try:
                os.mkdir(dirpath)
            except FileExistsError:
                pass
            pyfits.writeto(os.path.join(dirpath,'frequency4D.fits'),freq4D)

        if buffer_dec is not None:
            # Wait for buffer to end
            buffer_thread.join()
        
        return tn


    def acq_frequencies(self, fvec, nframes:int = 500):

        max_freq4D = self.interf.getFrameRate()
        fcy= int(max_freq4D/max(fvec))
        for freqPI in fvec:
            freq4D = freqPI*fcy  # same temporal sampling for each PI freq
            #freq = freqPI*nframes/Ncycles
            #max_freq4D = self.interf.getFrameRate()
            #freq4D = np.min(max_freq4D,int(freq*5)/5) # multiple of 0.2 [Hz]
            print(f'Piezo freq: {freqPI:1.0f} [Hz], 4D freq: {freq4D:1.1f} [Hz]')
            self.acq_freq(freqPI, freq4D, nframes = nframes)


    def acq_sync_sweep(self, fmin = 30, fmax = 110, duration = 11, ampPI = 2, 
                  nframes:int = 2250, produce:bool = False):

        # Generate new tn
        tn = self._generate_tn()

        # Start piezo and wait for steady state
        self.wavegen.set_wave(ch=1,ampl=ampPI,offs=0,freq=fmin,wave_form='SIN')
        time.sleep(4) 
        
        # Prepare for sweep
        self.wavegen.set_sweep(1,fmin,fmax,duration,amp=ampPI)
        self.interf.setTriggerMode(True)
        self.wavegen.set_burst(ch=2,freq=225,Ncycles=nframes)

        # Use multithreading for synchronization
        def trigger():
            self.wavegen.trigger_sweep(1)
        def capture():
            self.interf.capture(nframes,tn)

        wg_thread = Thread(target = trigger)
        interf_thread = Thread(target = capture)

        interf_thread.start() # run capture
        wg_thread.start() # trigger acquisition

        interf_thread.join() # wait for capture to end
        wg_thread.join() # wait for trigger to end (should be immediate)

        # Save and log data
        self.webdaq.stop_schedule_when_job_ends(job_id = 0) 
        self.interf.setTriggerMode(False)

        self.sync_and_save_last_wdfile(tn)

        if self.mx is not None:
            self.mx.save_pressure(self.config_path,tn)

        if self.dm is not None:
            self.dm.save_state(self.config_path,tn)

        if produce:
            self.interf.produce(tn)

        return tn


    @staticmethod
    def _generate_tn():
        tn = datetime.now()
        tnout = tn.strftime('%Y%m%d_%H%M%S')

        return tnout
    
    @staticmethod
    def sync_and_save_last_wdfile(tn):
        sp.wdsync()
        wdf = sp.last_wdfile()
        sp.savefile(wdf, tn)




