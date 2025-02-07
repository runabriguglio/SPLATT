from splattsw.devices.devices_scpi import SCPI


class PowerSupplier(SCPI):

    def __init__(self, IP:str='192.168.0.101', port:int=5555, TIMEOUT=2):

        super().__init__(IP,port,TIMEOUT)
        self.default_state = 'splatt.RSF'

    def reset(self):
        self.connect()
        self.rst()
        self.close()

    def switch_on(self, ch):
        self.connect()
        self.tx_txt(':OUTP CH'+str(ch)+',ON')
        self.close()

    def switch_off(self, ch):
        self.connect()
        self.tx_txt(':OUTP CH'+str(ch)+',OFF')
        self.close()

    def switch_state(self, ch):
        self.connect()
        state = self.txrx_txt(':OUTP? CH'+str(ch))
        state = 'ON' if state == '1' else 'OFF'
        self.close()
        return state

    def read_voltage(self, ch):
        self.connect()
        voltage = self.txrx_txt(':MEAS:VOLT? CH'+str(ch))
        self.close()
        return float(voltage)

    def read_current(self, ch):
        self.connect()
        current = self.txrx_txt(':MEAS:CURR? CH'+str(ch))
        self.close()
        return float(current)

    def load_default_state(self):
        self.connect()
        self.tx_txt(':MEM:LOAD '+ self.default_state)
        self.close()

    def load_saved_state(self, state_name):
        self.connect()
        self.tx_txt(':MEM:LOAD '+ state_name)
        self.close()

    def save_current_state(self, state_name):
        self.connect()
        self.tx_txt(':MEM:STOR '+ state_name + '.RSF')
        self.close()
