import matlab.engine
import Pyro4
import numpy as np

@Pyro4.expose
class MatlabEngine(object):

    def connect_matlab(self, shared_session_name: str = 'splattEngine'):
        names = matlab.engine.find_matlab()
        if shared_session_name in names:
            print(f"Connecting to '{shared_session_name}' ...")
            self.eng = matlab.engine.connect_matlab(shared_session_name)
        else:
            print(f"Shared MATLAB session '{shared_session_name}' not found, starting a new session ...")
            self.eng = matlab.engine.start_matlab()
            self.eng.cd(r'/home/labot/SPLATT_SW/Matlab_2024/Matlab/Scripts')
            self.eng.eval('splattInit',nargout=0)

    def send_command(self, cmd_str, wait4reply:bool = True):
        if wait4reply is False:
            self._oneway_command(cmd_str)
        else:
            self._command(cmd_str)

    def read_data(self, command_to_read_data:str, n_args_out: int = 1):
        # Pyro does not seem to support numpy: convert any arrays after the call
        mat_data = self.eng.eval(str(command_to_read_data),nargout=n_args_out)

        data = mat_data
        if n_args_out > 1:
            data = []
            for val in mat_data:
                data.append(val)

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
