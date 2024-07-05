#!/usr/bin/python
# python gen_wave.py IP ampl offset freq wave_form
 
import sys
import time
import  SPLATT.splattsw.devices.devices_scpi as scpi

print("Warning!! Load your device with function ''update_device(deviceName)''")

global name

def update_device(deviceName):
    global name
    name = deviceName
    
ampTrigger4d = 2
offsetTrigger4d = 2.3
ampTrigger4dPulse = 4.
offsetTrigger4dPulse = 5
pulseDurationTrigger = 100e-6
amplifierGain = 10.
ampPiezo = 2.

def splatt_trigger( freqPI, freq4d, ampPI=ampPiezo):
    if ampPI > 4:
        raise Exception('Amplitude too large!!')

    set_wave1(ampPI/amplifierGain, 0, freqPI, 'SINE')
    trigg4d(freq4d)
    phase_align()
    
def splatt_trigger2( freqPI, freq4d, ampPI=ampPiezo):
    #if ampPI > 4:
    #    raise Exception('Amplitude too large!!')

    set_wave1(ampPI/amplifierGain, 0, freqPI, 'SINE')
    trigg4d_pulse(freq4d)
    phase_align()

   
def splatt_pulse(freqPI, freq4d, ampPI=ampPiezo):
    if ampPI > 4:
        raise Exception('Amplitude too large!!')
    trigg4d_pulse(freq4d)
    pulse_train(1, ampPI/amplifierGain,0,freqPI,0.05)

 
def trigg4d_pulse(freq):
    amp=ampTrigger4dPulse/amplifierGain
    offs = offsetTrigger4dPulse/amplifierGain
    dc = pulseDurationTrigger*freq
    set_pulse(2, amp, offs, freq, dc)

def trigg4d(freq):
    set_wave2(ampTrigger4d/amplifierGain, offsetTrigger4d/amplifierGain, freq, 'SQUARE')    

def single_pulse(ch):
    pulseamp = 2
    '''
    time.sleep(0.2)
    clear_wave(ch)
    time.sleep(0.2)
    wave_on(ch)
    '''
    set_pulse(ch, pulseamp/amplifierGain, pulseamp/amplifierGain, 20, 0.1)
    time.sleep(0.2)
    clear_wave(ch)
    wave_on(2)


def pulse_train(ch, amp, offs, freq, dc, duration=None):
    clear_wave(ch)
    time.sleep(0.2)
    wave_on(ch)
    set_pulse(ch, amp, offs, freq, dc)
    print('waiting...')
    if duration != None:
        
        time.sleep(duration)
        clear_wave(ch)    
    else:
        print('No duration set, waiting for clear_wave')

def set_waves(c0,c1,f0,f1,ph, type0='SINE', type1='SINE'):
    set_wave1(c0,0,f0,type0)
    set_wave2(c1,0,f1,type1)



def blink():
    rp_s = scpi.scpi(name)

    if (len(sys.argv) > 2):
        led = int(sys.argv[2])
    else:
        led = 0

    print ("Blinking LED["+str(led)+"]")

    period = 1 # seconds

    for i in range(60):
        time.sleep(period/2.0)
        rp_s.tx_txt('DIG:PIN LED' + str(led) + ',' + str(1))
        time.sleep(period/2.0)
        rp_s.tx_txt('DIG:PIN LED' + str(led) + ',' + str(0))

def set_wave1(ampl, offs, freq, wave_form):
    rp_s = scpi.scpi(name)

    rp_s.tx_txt('SOUR1:FUNC ' + str(wave_form).upper())
    rp_s.tx_txt('SOUR1:FREQ:FIX ' + str(freq))
    rp_s.tx_txt('SOUR1:VOLT ' + str(ampl))
    rp_s.tx_txt('SOUR1:VOLT:OFFS ' + str(offs))

    #Enable output
    rp_s.tx_txt('OUTPUT1:STATE ON')
    rp_s.close()
    
def set_wave2(ampl, offs, freq, wave_form):
    rp_s = scpi.scpi(name)

    rp_s.tx_txt('SOUR2:FUNC ' + str(wave_form).upper())
    rp_s.tx_txt('SOUR2:FREQ:FIX ' + str(freq))
    rp_s.tx_txt('SOUR2:VOLT ' + str(ampl))
    rp_s.tx_txt('SOUR2:VOLT:OFFS ' + str(offs))

    #Enable output
    rp_s.tx_txt('OUTPUT2:STATE ON')
    rp_s.close()
 
def clear_wave(ch):
    rp_s = scpi.scpi(name)
    rp_s.tx_txt('OUTPUT'+str(ch)+':STATE OFF')
    rp_s.tx_txt('SOUR'+str(ch)+':FREQ:FIX 0')
    rp_s.tx_txt('SOUR'+str(ch)+':VOLT 0')
    rp_s.tx_txt('SOUR'+str(ch)+':VOLT:OFFS 0')
    rp_s.close()


def clear_wave1():
    rp_s = scpi.scpi(name)
    #Enable output
    rp_s.tx_txt('OUTPUT1:STATE OFF')
    rp_s.tx_txt('SOUR1:FREQ:FIX 0')
    rp_s.tx_txt('SOUR1:VOLT 0')
    rp_s.tx_txt('SOUR1:VOLT:OFFS 0')
    rp_s.close()


def clear_wave2():
    rp_s = scpi.scpi(name)
    #Enable output
    rp_s.tx_txt('OUTPUT2:STATE OFF')
    rp_s.tx_txt('SOUR2:FREQ:FIX 0')
    rp_s.tx_txt('SOUR2:VOLT 0')
    rp_s.tx_txt('SOUR2:VOLT:OFFS 0')
    #rp_s.tx_txt('OUTPUT2:STATE OFF')
    rp_s.close()

def wave_on(ch):
    rp_s = scpi.scpi(name)
    rp_s.tx_txt('OUTPUT'+str(ch)+':STATE ON')
    rp_s.close()

def wave_off(ch):
    rp_s = scpi.scpi(name)
    rp_s.tx_txt('OUTPUT'+str(ch)+':STATE OFF')
    rp_s.close()

def waves_off():
    rp_s = scpi.scpi(name)
    #Enable output
    rp_s.tx_txt('OUTPUT:STATE OFF')
    rp_s.close()

def waves_on():
    rp_s = scpi.scpi(name)
    #Enable output
    rp_s.tx_txt('OUTPUT:STATE ON')
    rp_s.close()

def generator_reset():#reset to default settings: 5V!!!! 
    rp_s = scpi.scpi(name)
    rp_s.tx_txt('GEN:RST')
    rp_s.close()

def phase_align():
    time.sleep(0.1)
    rp_s = scpi.scpi(name)
    rp_s.tx_txt('PHAS:ALIGN')
    rp_s.close()
    time.sleep(0.1)

def set_phase(ch, phase):
    phase_align()
    rp_s = scpi.scpi(name)
    sch = str(ch)
    ss = str(phase)
    scomm = 'SOUR'+sch+':PHAS '+ss
    rp_s.tx_txt(scomm)

    rp_s.close()

def set_pulse(ch, amp, offs, freq, dc):
    rp_s = scpi.scpi(name)
    sch = str(ch)
    samp = str(amp)
    soffs = str(offs)
    sfreq = str(freq)
    sdc = str(dc)
    rp_s.tx_txt('SOUR'+sch+':FUNC PWM')
    rp_s.tx_txt('SOUR'+str(ch)+':FREQ:FIX ' + str(freq))
    
    
    rp_s.tx_txt('SOUR'+sch+':VOLT ' + str(amp))
    rp_s.tx_txt('SOUR'+sch+':VOLT:OFFS ' + str(offs))
    rp_s.tx_txt('SOUR'+sch+':DCYC ' + str(dc))
    rp_s.close()



if __name__=="__main__":
    wave_form = 'sine'
    freq = 20
    ampl = 0.1
    offs=0.1
    setwave1(ampl, offs, freq, wave_form)
