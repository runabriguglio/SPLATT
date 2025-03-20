import matlab.engine
import Pyro4
import numpy as np

@Pyro4.expose
class MatlabEngine(object):

    def start_engine(self):
        self.eng = matlab.engine.connect_matlab()
        self.eng.cd(r'/home/labot/SPLATT_SW/Matlab_2024/Matlab/Scripts')

    def send_command(self, cmd_str, wait4reply:bool = True):
        if wait4reply is False:
            self._oneway_command(self,cmd_str)
        else:
            self._oneway_command(self,cmd_str)

    def get_data(self, command_to_read_data:str, n_args_out: int = 1):
        # Note that Pyro does not seem to support numpy, convert any arrays after the call
        mat_data = self.eng.eval(str(command_to_read_data),nargout=n_args_out)
        data = []
        for n in range(n_args_out):
            val = np.array(mat_data[n])
            data.append(val.tolist())
        return data

    def stop_engine(self):
        self.eng.quit()

    def _command(self,command):
        self.eng.eval(str(command)+';', nargout=0)

    @Pyro4.oneway
    def _oneway_command(self,command):
        self.eng.eval(str(command)+';', nargout=0)


def main():

    matlab_eng = MatlabEngine()

    Pyro4.Daemon.serveSimple( {matlab_eng: "matlab_engine"},
             host="193.206.155.220", port=9090, ns=False, verbose=True)


if __name__=="__main__":
    main()
