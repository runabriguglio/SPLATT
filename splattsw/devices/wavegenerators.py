#!/usr/bin/python
# python gen_wave.py IP ampl offset freq wave_form

import sys
import time
from splattsw.devices.devices_scpi import SCPI

class WaveGenerator(SCPI):

    def __init__(self, IP:str='192.168.0.100', port:int=5555, TIMEOUT=2):

        super().__init__(IP,port,TIMEOUT)
        # self.name = device_name
        # self = scpi.scpi(self.name)
        self.ampTrigger4d = 2
        self.offsetTrigger4d = 2.3
        self.ampTrigger4dPulse = 4.
        self.offsetTrigger4dPulse = 5
        self.pulseDurationTrigger = 100e-6
        self.amplifierGain = 10.
        self.ampPiezo = 2.
        self.delay = 0.1


    def splatt_trigger(self, freqPI, freq4d, ampPI=None):
        if ampPI > 4:
            raise Exception('Amplitude too large!!')
        elif ampPI is None:
            ampPI = self.ampPiezo

        self.set_wave1(ampPI/self.amplifierGain, 0, freqPI, 'SIN')
        self.trigg4d(freq4d)
        self.phase_align()

    def splatt_trigger2(self, freqPI, freq4d, ampPI=None):
        #if ampPI > 4:
        #    raise Exception('Amplitude too large!!')
        if ampPI is None:
            ampPI = self.ampPiezo
        self.set_wave1(ampPI/self.amplifierGain, 0, freqPI, 'SIN')
        self.trigg4d_pulse(freq4d)
        self.phase_align()


    def splatt_pulse(self, freqPI, freq4d, ampPI=None):
        if ampPI > 4:
            raise Exception('Amplitude too large!!')
        elif ampPI is None:
            ampPI = self.ampPiezo
        self.trigg4d_pulse(freq4d)
        self.pulse_train(1, ampPI/self.amplifierGain,0,freqPI,0.05)


    def trigg4d_pulse(self, freq):
        amp=self.ampTrigger4dPulse/self.amplifierGain
        offs = self.offsetTrigger4dPulse/self.amplifierGain
        dc = self.pulseDurationTrigger*freq
        self.set_pulse(2, amp, offs, freq, dc)

    def trigg4d(self, freq):
        self.set_wave2(
            self.ampTrigger4d/self.amplifierGain,
            self.offsetTrigger4d/self.amplifierGain,
            freq,
            'SQUARE'
        )

    # def reset(self):
    #     self.connect()
    #     print("WARNING! Following commmands might be corrupted")
    #     self.rst()
    #     self.close()

    def state_on(self, ch):
        self.wave_on(ch)

    def single_pulse(self, ch):
        pulseamp = 2
        '''
        time.sleep(0.2)
        clear_wave(ch)
        time.sleep(0.2)
        wave_on(ch)
        '''
        self.set_pulse(ch, pulseamp/self.amplifierGain, pulseamp/self.amplifierGain, 20, 0.1)
        time.sleep(0.2)
        self.clear_wave(ch)
        self.wave_on(ch)

    def sweep(self, ch, fmin, fmax, period=20, amp=1):
        self.connect()
        if fmin<10 or fmax>150:
            raise ValueError(f"Are you sure you want to go from {fmin} Hz to {fmax} Hz ?!")
        self.tx_txt(f"OUTPUT{ch}:STATE ON")
        self.tx_txt(f":SOUR{ch}:FREQ:STAR {fmin}")
        self.tx_txt(f":SOUR{ch}:FREQ:STOP {fmax}")
        self.tx_txt(f":SOUR{ch}:SWE:TIME {period}")
        self.tx_txt(f":SOUR{ch}:SWE:RTIM 0")
        self.tx_txt(f":SOUR{ch}:VOLT {amp}") # Set amplitude
        self.tx_txt(':SOUR{ch}:SWE:TRIG:SOUR MAN')
        #self.tx_txt(f":SOUR{ch}:SWE:TRIG") # Start sweep immediately
        self.close()

    def pulse_train(self, ch, amp, offs, freq, dc, duration=None):
        self.clear_wave(ch)
        time.sleep(0.2)
        self.wave_on(ch)
        self.set_pulse(ch, amp, offs, freq, dc)
        print('waiting...')
        if duration != None:
            time.sleep(duration)
            self.clear_wave(ch)
        else:
            print('No duration set, waiting for clear_wave')

    def set_waves(self, c0,c1,f0,f1, type0='SIN', type1='SIN'):
        self.set_wave1(c0,0,f0,type0)
        self.set_wave2(c1,0,f1,type1)


    def blink(self):
        self.connect()
        if (len(sys.argv) > 2):
            led = int(sys.argv[2])
        else:
            led = 0
        print ("Blinking LED["+str(led)+"]")
        period = 1 # seconds
        for i in range(60):
            time.sleep(period/2.0)
            self.tx_txt('DIG:PIN LED' + str(led) + ',' + str(1))
            time.sleep(period/2.0)
            self.tx_txt('DIG:PIN LED' + str(led) + ',' + str(0))

    def set_wave(self, ch, ampl, offs, freq, wave_form):
        self.connect()
        self.tx_txt('SOUR'+str(ch)+':FUNC ' + str(wave_form).upper())
        time.sleep(self.delay)
        self.tx_txt('SOUR'+str(ch)+':FREQ:FIX ' + str(freq))
        time.sleep(self.delay)
        self.tx_txt('SOUR'+str(ch)+':VOLT ' + str(ampl))
        time.sleep(self.delay)
        self.tx_txt('SOUR'+str(ch)+':VOLT:OFFS ' + str(offs))
        time.sleep(self.delay)
        #Enable output
        self.tx_txt('OUTPUT'+str(ch)+':STATE ON')
        self.close()


    def set_wave1(self, ampl, offs, freq, wave_form):
        self.set_wave(1, ampl, offs, freq, wave_form)

    def set_wave2(self, ampl, offs, freq, wave_form):
        self.set_wave(2, ampl, offs, freq, wave_form)

    def clear_wave(self, ch):
        self.connect()
        self.tx_txt(f'OUTPUT{ch}:STATE OFF')
        self.tx_txt(f'SOUR{ch}:VOLT:OFFS 0')
        self.tx_txt(f":SOUR{ch}:PHAS 0")
        self.tx_txt(f":SOUR{ch}:SWE:STAT 0") # Sweep mode off
        self.close()


    def clear_wave1(self):
        self.clear_wave(1)


    def clear_wave2(self):
        self.clear_wave(2)

    def wave_on(self, ch):
        self.connect()
        self.tx_txt('OUTPUT'+str(ch)+':STATE ON')
        self.close()

    def wave_off(self, ch):
        self.connect()
        self.tx_txt('OUTPUT'+str(ch)+':STATE OFF')
        self.close()

    def waves_off(self):
        self.connect()
        #Enable output
        self.tx_txt('OUTPUT:STATE OFF')
        self.close()

    def waves_on(self):
        self.connect()
        #Enable output
        self.tx_txt('OUTPUT:STATE ON')
        self.close()

    def phase_align(self):
        time.sleep(0.1)
        self.connect()
        self.tx_txt('PHAS:ALIGN')
        self.close()
        time.sleep(0.1)

    def set_phase(self, ch, phase):
        self.phase_align()
        self.connect()
        sch = str(ch)
        ss = str(phase)
        scomm = 'SOUR'+sch+':PHAS '+ss
        self.tx_txt(scomm)
        self.close()

    def set_pulse(self, ch, amp, offs, freq, dc):
        self.connect()
        sch = str(ch)
        samp = str(amp)
        soffs = str(offs)
        sfreq = str(freq)
        sdc = str(dc)
        self.tx_txt('SOUR'+sch+':FUNC PWM'); time.sleep(self.delay)
        self.tx_txt('SOUR'+ sch +':FREQ:FIX ' + sfreq); time.sleep(self.delay)
        self.tx_txt('SOUR'+sch+':VOLT ' + samp); time.sleep(self.delay)
        self.tx_txt('SOUR'+sch+':VOLT:OFFS ' + soffs); time.sleep(self.delay)
        self.tx_txt('SOUR'+sch+':DCYC ' + sdc); time.sleep(self.delay)
        self.close()



if __name__=="__main__":
    wg = WaveGenerator()
    wave_form = 'sin'
    freq = 20
    ampl = 0.1
    offs= 0
    wg.setwave1(ampl, offs, freq, wave_form)
