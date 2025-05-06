
#from splattsw import splatt_analysis as sp
from splattsw import acceleration_analysis as sp
from splattsw.devices.wavegenerators import WaveGenerator as wg
from splattsw.devices.webDAQ import WebDAQ as wdq

#import aoptics
#pyconf = '/mnt/libero/SPLATTData/Data/SysConfigurations/configuration.yaml'
#aoptics.load_configuration_file(pyconf)

from aoptics.devices.interferometer import PhaseCam
from splattsw.devices.moxa_io import Moxa_ai0
from splattsw.devices.deformable_mirror import SPLATTEngine

import time
import numpy as np
from datetime import datetime

class Acquisition():

    def __init__(self):

        # Define interferometer
        self.interf = None
        try:
            self.interf=PhaseCam('4020')
        except:
            print('Interferometer not found')

        # Accelerometers
        self.webdaq=wdq()
        self.webdaq.connect()

        # Signal generator
        self.wavegen = wg()

        # SPLATT dm
        self.dm = None
        try:
            dm = SPLATTEngine()
            self.dm = dm._eng
        except:
            print('SPLATT dm not found')

        # MOXA
        self.mx = None
        try:
            self.mx = Moxa_ai0()
        except:
            print('Moxa device not found')



    def acq_sweep(self, fmin = 30, fmax = 110, duration = 11, ampPI = 2, nframes:int = 2250, chPI:int = 1, produce:bool = False):

        # Generate new tn
        tn = self._generate_tn()

        # Setup sweep parameters
        self.wavegen.set_wave(ch=chPI,ampl=ampPI,offs=0,freq=fmin,wave_form='SIN')
        time.sleep(4) # wait for steady state
        self.wavegen.set_sweep(chPI,fmin,fmax,duration,amp=ampPI)
        #self.wavegen.set_burst(ch=2,freq=225,Ncycles=nframes)

        # Start acquisition and sweep
        self.webdaq.start_schedule()
        self.wavegen.trigger_sweep(chPI)
        time.sleep(0.5)
        self.interf.capture(nframes, tn)

        # Wait for webDAQ acquisition to end
        self.webdaq.stop_schedule_when_job_ends(job_id = 0)

        # Post-processing
        sp.wdsync()
        wdfile=sp.last_wdfile()
        sp.savefile(wdfile,tn) # saving to fits while renaming with correct tn

        if self.mx is not None:
            pres = self.mx.read_pressure()
            print(f'The vacuarium pressure is {pres:1.3f} [bar]')

        if produce:
            self.interf.produce(tn)

        return tn



    def acq_freq(self, freqPI,  ampPI=2, nframes:int = 2000, chPI:int = 1, produce:bool = False, buffer:bool = False):

        # Generate new tn
        tn = self._generate_tn()

        # Start acquisition and sine wave
        self.wavegen.set_wave(ch=chPI, ampl=ampPI, offs=0, freq=freqPI, wave_form='SIN')
        if buffer:
            self.dm.send(f'clear opts; opts.dec = 0; opts.sampleNR = 256; opts.save2mat = 0; opts.save2fits = 1; opts.tn = {tn}')
            self.eng.oneway_send("splattAcqBufInt({'sabi32_Distance','sabi32_pidCoilOut'},opts)")
        t0 = time.time()
        time.sleep(4) # wait for transient
        self.webdaq.start_schedule()
        if nframes > 0:
            self.interf.capture(nframes, tn)
        
        # Wait for webDAQ acquisition to end
        self.webdaq.stop_schedule_when_job_ends(job_id = 0)

        # Acquisition end
        self.wavegen.wave_off(chPI)
        
        # Post-processing
        time.sleep(1)
        sp.wdsync()
        wdfile=sp.last_wdfile()
        sp.savefile(wdfile,tn) # saving to fits while renaming with correct tn

        if self.mx is not None:
            p_bar = self.mx.read_pressure()
            print(f'The vacuarium pressure is {p_bar:1.3f} [bar]')

        if np.logical_and( nframes > 0, produce):
            self.interf.produce(tn)

        if buffer:
            # Wait for buffer to end before reading tn
            t1 = time.time()
            dt = t1-t0

            if dt < 45:
                time.sleep(45-dt)
        
        return tn


    @staticmethod
    def _generate_tn():
        tn = datetime.now()
        tn.strftime('%Y%m%d_%H%M%S')

        return tn




