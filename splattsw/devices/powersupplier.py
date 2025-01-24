import splattsw.devices.devices_scpi as scpi


class PowerSupplier:

    def __init__(self, device_name:str='Rigol_PowerSupplier'):
        self.name = device_name
        self.TIMEOUT = 2
        self.rp_s = scpi.scpi(self.name, self.TIMEOUT)
        self.default_state = 'splatt.RSF'

    def reset(self):
        self.rp_s.connect()
        self.rp_s.rst()
        self.rp_s.close()

    def switch_on(self, ch):
        self.rp_s.connect()
        self.rp_s.tx_txt(':OUTP CH'+str(ch)+',ON')
        self.rp_s.close()

    def switch_off(self, ch):
        self.rp_s.connect()
        self.rp_s.tx_txt(':OUTP CH'+str(ch)+',OFF')
        self.rp_s.close()

    def switch_state(self, ch):
        self.rp_s.connect()
        state = self.rp_s.txrx_txt(':OUTP? CH'+str(ch))
        state = 'ON' if state == '1' else 'OFF'
        self.rp_s.close()
        return state

    def read_voltage(self, ch):
        self.rp_s.connect()
        voltage = self.rp_s.txrx_txt(':MEAS:VOLT? CH'+str(ch))
        self.rp_s.close()
        return float(voltage)

    def read_current(self, ch):
        self.rp_s.connect()
        current = self.rp_s.txrx_txt(':MEAS:CURR? CH'+str(ch))
        self.rp_s.close()
        return float(current)

    def load_saved_state(self, state_name):
        self.rp_s.connect()
        if state_name is None:
            state_name = self.default_state
        self.rp_s.tx_txt(':MEM:LOAD '+ state_name)
        self.rp_s.close()

    def save_current_state(self, state_name):
        self.rp_s.connect()
        self.rp_s.tx_txt(':MEM:STOR '+ state_name + '.RSF')
        self.rp_s.close()
