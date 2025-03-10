import matlab.engine
import Pyro4
#import numpy as np

@Pyro4.expose
class MatlabEngine(object):

    def start_engine(self):
        self.eng = matlab.engine.start_matlab()
        #self.eng.desktop(nargout=0)
        self.eng.cd(r'/home/labot/SPLATT_SW/Matlab_2024/Matlab/Scripts')

    def send_command(self,command):
        self.eng.eval(str(command)+';', nargout=0)

    def get_data(self, command_to_read_data:str):
        # Note that Pyro does not seem to support numpy, convert
        # any arrays after the call
        data = self.eng.eval(str(command_to_read_data),nargout=1)
        return data

    def stop_engine(self):
        self.eng.quit()


def main():

    matlab_eng = MatlabEngine()

    Pyro4.Daemon.serveSimple( {matlab_eng: "matlab_engine"},
             host="193.206.155.220", port=9090, ns=False, verbose=True)


if __name__=="__main__":
    main()
