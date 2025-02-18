from devices.moxa_class import Moxa

# Inherited classes: add your own moxa!
class Moxa_pt0(Moxa):

    def __init__(self, ip = '193.206.155.40', nchannels:int = 6,
                api_addr_ext = '/api/slot/0/io/rtd',
                valuestr = 'rtdValueScaled', valueid = 'rtd'):

        super().__init__(ip, nchannels, api_addr_ext, valuestr, valueid)



class Moxa_pt1(Moxa):

    def __init__(self, ip = '193.206.155.41', nchannels:int = 6,
                api_addr_ext = '/api/slot/0/io/rtd',
                valuestr = 'rtdValueScaled', valueid = 'rtd'):

        super().__init__(ip, nchannels, api_addr_ext, valuestr, valueid)



class Moxa_ai0(Moxa):

    def __init__(self, ip = '193.206.155.47', nchannels:int = 8,
                api_addr_ext = '/api/slot/0/io/ai',
                valuestr = 'aiValueScaled', valueid = 'ai'):

        super().__init__(ip, nchannels, api_addr_ext, valuestr, valueid)

    def read_pressure(self):
        data = self.read()
        pres = data[6]
        return pres



class Moxa_di0(Moxa):

    def __init__(self, ip = '193.206.155.141', nchannels:int = 16,
                api_addr_ext = '/api/slot/0/io/di',
                valuestr = 'diValueScaled', valueid = 'di'):

        super().__init__(ip, nchannels, api_addr_ext, valuestr, valueid)