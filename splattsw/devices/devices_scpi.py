"""SCPI access to Red Pitaya or Rigol devices."""
import time
import socket


class SCPI(object):
    """SCPI class used to access Red Pitaya over an IP network."""
    delimiter = '\r\n'

    def __init__(self, IP, port, timeout = None):
        """Initialize object and open IP connection.
        Host IP should be a string in parentheses, like '192.168.1.100'.
        """
        # self.name = deviceName
        self.cip = IP
        self.port = port
        self.timeout = timeout
        # self._select_device()
        self.connect()
        # try:
        #     self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        #     if self.timeout is not None:
        #         self._socket.settimeout(self.timeout)

        #     self._socket.connect((self.cip, self.port))

        # except socket.error as e:
        #     print('SCPI >> connect({:s}:{:d}) failed: {:s}'.format(self.cip, self.port, e))

    def __del__(self):
        if self._socket is not None:
            self._socket.close()
        self._socket = None

    def connect(self):
        """Connect to device IP."""
        try:
            self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            if self.timeout is not None:
                self._socket.settimeout(self.timeout)
            self._socket.connect((self.cip, self.port))
        except socket.error as e:
            print(f"SCPI >> connect({self.cip}:{self.port}) failed: '{e}'")

    def close(self):
        """Close IP connection."""
        self.__del__()

    def rx_txt(self, chunksize = 4096):
        """Receive text string and return it after removing the delimiter."""
        msg = ''
        while 1:
            chunk = self._socket.recv(chunksize + len(self.delimiter)).decode('utf-8') # Receive chunk size of 2^n preferably
            msg += chunk
            if (len(chunk) and chunk[-2:] == self.delimiter):
                break
        return msg[:-2]

    def rx_arb(self):
        numOfBytes = 0
        """ Recieve binary data from scpi server"""
        str=b''
        while (len(str) != 1):
            str = (self._socket.recv(1))
        if not (str == b'#'):
            return False
        str=b''
        while (len(str) != 1):
            str = (self._socket.recv(1))
        numOfNumBytes = int(str)
        if not (numOfNumBytes > 0):
            return False
        str=b''
        while (len(str) != numOfNumBytes):
            str += (self._socket.recv(1))
        numOfBytes = int(str)
        str=b''
        while (len(str) != numOfBytes):
            str += (self._socket.recv(4096))
        return str

    def tx_txt(self, msg):
        """Send text string ending and append delimiter."""
        return self._socket.send((msg + self.delimiter).encode('utf-8'))

    # def txrx_txt(self, msg):
    #     """Send/receive text string."""
    #     self.tx_txt(msg)
    #     return self.rx_txt()

    def txrx_txt(self, msg:str, cs:int=4096):
        """Send a reading message and outputs the response and the error message."""
        self._socket.send((msg + self.delimiter).encode('utf-8'))
        time.sleep(0.2)
        response = self._socket.recv(cs).decode('utf-8').strip()
        self._socket.send((":SYST:ERR?" + self.delimiter).encode('utf-8'))
        time.sleep(0.2)
        error_response = self._socket.recv(cs).decode('utf-8').strip()
        if error_response != '0,"No error"':
            print(f"ERROR: {error_response}")
        return response



# IEEE Mandated Commands

    def cls(self):
        """Clear Status Command"""
        return self.tx_txt('*CLS')

    def ese(self, value: int):
        """Standard Event Status Enable Command"""
        return self.tx_txt('*ESE {}'.format(value))

    def ese_q(self):
        """Standard Event Status Enable Query"""
        return self.txrx_txt('*ESE?')

    def esr_q(self):
        """Standard Event Status Register Query"""
        return self.txrx_txt('*ESR?')

    def idn_q(self):
        """Identification Query"""
        return self.txrx_txt('*IDN?')

    def opc(self):
        """Operation Complete Command"""
        return self.tx_txt('*OPC')

    def opc_q(self):
        """Operation Complete Query"""
        return self.txrx_txt('*OPC?')

    def rst(self):
        """Reset Command"""
        return self.tx_txt('*RST')

    def sre(self):
        """Service Request Enable Command"""
        return self.tx_txt('*SRE')

    def sre_q(self):
        """Service Request Enable Query"""
        return self.txrx_txt('*SRE?')

    def stb_q(self):
        """Read Status Byte Query"""
        return self.txrx_txt('*STB?')

# :SYSTem

    def err_c(self):
        """Error count."""
        return self.txrx_txt('SYST:ERR:COUN?')

    def err_c(self):
        """Error next."""
        return self.txrx_txt('SYST:ERR:NEXT?')
