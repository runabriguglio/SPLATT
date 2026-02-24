#from splattsw.devices import moxa_io as mx

class Temperature:
    """
    HOW TO USE IT::

    from m4.configuration import start
    ott, interf = start.create_ott(conf='.../youConf.yaml')
    from m4.mini_OTT.measurements import Measurements
    meas = Measurements(ott, interf)
    """

    def __init__(self, datalogger):
        """The constructor"""
        self._datalogger = datalogger


    def getTemperature():
        tempv = datalogger.read()
        return tempv


