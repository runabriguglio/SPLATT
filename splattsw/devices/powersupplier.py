import splattsw.devices.devices_scpi as scpi

name = 'Rigol_PowerSupplier'
TIMEOUT = 2
default_state = 'splatt.RSF'
# USE THE FOLLOWING IF MORE THAN 1 POWER SUPPLIERS ARE AVAILABLE
# global name
# print("Warning!! Load your device with function ''update_device(deviceName)''")
# def update_device(deviceName):
#     global name
#     name = deviceName

def reset():
    rp_s = scpi.scpi(name, TIMEOUT)
    rp_s.rst()
    rp_s.close()

def switch_on(ch):
    rp_s = scpi.scpi(name, TIMEOUT)
    rp_s.tx_txt(':OUTP CH'+str(ch)+',ON')
    rp_s.close()

def switch_off(ch):
    rp_s = scpi.scpi(name, TIMEOUT)
    rp_s.tx_txt(':OUTP CH'+str(ch)+',OFF')
    rp_s.close()

def switch_state(ch):
    rp_s = scpi.scpi(name, TIMEOUT)
    state = rp_s.txrx_txt(':OUTP? CH'+str(ch))
    state = 'ON' if state == '1' else 'OFF'
    rp_s.close()
    return state

def read_voltage(ch):
    rp_s = scpi.scpi(name, TIMEOUT)
    voltage = rp_s.txrx_txt(':MEAS:VOLT? CH'+str(ch))
    rp_s.close()
    return float(voltage)

def read_current(ch):
    rp_s = scpi.scpi(name, TIMEOUT)
    current = rp_s.txrx_txt(':MEAS:CURR? CH'+str(ch))
    rp_s.close()
    return float(current)

def load_saved_state(state_name = default_state):
    rp_s = scpi.scpi(name, TIMEOUT)
    rp_s.tx_txt(':MEM:LOAD '+ state_name)
    rp_s.close()

def save_current_state(state_name):
    rp_s = scpi.scpi(name, TIMEOUT)
    rp_s.tx_txt(':MEM:STOR '+ state_name + '.RSF')
    rp_s.close()
