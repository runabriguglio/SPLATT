import SPLATT.splattsw.devices.devices_scipi as scpi

print("Warning!! Load your device with function ''update_device(deviceName)''")

global name

def update_device(deviceName):
    global name
    name = deviceName

default_state = 'splatt.RSF'

def reset():
    print("WARNING! Following commmands might be corrupted")
    rp_s = scpi.scpi(name)
    rp_s.rst()
    rp_s.close()

def switch_on(ch):
    rp_s = scpi.scpi(name)
    rp_s.tx_txt(':OUTP CH'+str(ch)+' ON')
    rp_s.close()

def switch_off(ch):
    rp_s = scpi.scpi(name)
    rp_s.tx_txt(':OUTP CH'+str(ch)+' OFF')
    rp_s.close()

def load_saved_state(state_name):
    rp_s = scpi.scpi(name)
    rp_s.tx_txt(':MEM:LOAD '+ state_name)
    rp_s.close()