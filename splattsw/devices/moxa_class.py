# author: a5892731
# date: 01.06.2022
# last update: 14.06.2022
# version: 1.0

# description: This is REST API GET function for MOXA ioLogik E1212 device

'''
source:
https://github.com/zentec/moxa-iologic-1200-monitor/blob/master/monitor.py

'''

from requests import get
from json import loads
from time import sleep, time
import numpy as np


class Moxa():

    def __init__(self, ip, n_channels: int, api_addr_ext, value_str, value_id):

        self.ip = ip
        self.nchannels = n_channels
        self.api_address_extension = api_addr_ext
        self.valuestring = value_str
        self.valueId = value_id


    def getIP(self):
        return self.ip


    def read(self):
        address2use = 'http://'+self.ip
        out =get(address2use+self.api_address_extension,{"rel_rhy": "network"},headers={"Content-Type": "application/json", "Accept": "vdn.dac.v1"},timeout=3,)
        json_blob = loads(out.content.decode('utf-8'))
        vecout = []
        for i in range(0,self.nchannels):
            vecout.append(json_blob['io'][self.valueId][i][self.valuestring])
        vecout = np.array(vecout)
        return vecout


    def get_data(self, api_address_extension = "/api/slot/0/io/do"):
        try:
            self.reply = get(self.address + api_address_extension,
                             {"rel_rhy": "network"},
                             headers={"Content-Type": "application/json", "Accept": "vdn.dac.v1"},
                             timeout=0.1,
                            )

            self.replay_status_code = self.reply.status_code
            self.connection_error = False
            self.last_get_time = time()

        except:
            self.connection_error = True




'''
def read_ai():
    address2use = 'http://'+ip_moxaAI1240
    out =get(address2use+api_address_extension_ai,{"rel_rhy": "network"},headers={"Content-Type": "application/json", "Accept": "vdn.dac.v1"},timeout=0.1,)
    (json_blob['io']['ai'])[0]['aiBurnoutValue']

class GetRequestData():
    def __init__(self, address = "http://193.206.155.47"):
        self.address = address

        self.connection_error = None # True if there is no connection for more than 1 s
        self.data_error = None
        self.replay_status_code = None # connection statuses = 200; 404 etc
        self.reply = None # request response
        self.json_data = None # converted response to json
        """downloaded data >>>"""
        self.deviceUpTime = None

        self.DO = None
        self.DI = None

        self.last_get_time = time()

    def get_sysInfo(self):
        self.get_data(api_address_extension = "/api/slot/0/sysInfo/device")
        if not self.connection_error:
            self.convert_received_data()
            self.deviceUpTime = self.json_data['sysInfo']['device'][0]['deviceUpTime']

    def get_DI(self):
        self.get_data(api_address_extension="/api/slot/0/io/di")
        if not self.connection_error:
            self.convert_received_data()
            self.DI = self.json_data['io']['di']

    def get_DO(self):
        self.get_data(api_address_extension="/api/slot/0/io/do")
        if not self.connection_error:
            self.convert_received_data()
            self.DO = self.json_data['io']['do']

    def run(self):
        # get device status
        self.get_sysInfo()

        # get digital inputs
        self.get_DI()

        # get digital outputs
        self.get_DO()

    def convert_received_data(self):
        if not self.connection_error:
            try:
                json_blob = loads(self.reply.content.decode('utf-8'))
                self.json_data = json_blob

                self.data_error = False

            except ValueError:
                self.data_error = True



if __name__ == "__main__":
    data = GetRequestData()

    while True:
        data.run()

        sleep(1)
        print(data.deviceUpTime)
        print("di {}".format(data.DO))

        print("do {}".format(data.DO))

'''
