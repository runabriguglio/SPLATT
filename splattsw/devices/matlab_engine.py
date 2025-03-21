import matlab.engine
import Pyro4
import numpy as np

@Pyro4.expose
class MatlabEngine(object):

    def __init__(self, shared_session_name: str = 'splattEngine'):
        names = matlab.engine.find_matlab()

        if shared_session_name in names:
            print(f"Connecting to '{shared_session_name}' ...")
            self.eng = matlab.engine.connect_matlab(shared_session_name)
        else:
            print(f"No shared MATLAB session found, starting a new shared session '{shared_session_name}'...")
            self.eng = matlab.engine.connect_matlab()
            self.eng.cd(r'/home/labot/SPLATT_SW/Matlab_2024/Matlab/Scripts')
            self.eng.eval('splattInit',nargout=0)
            self.eng.shareEngine(f'{shared_session_name}') # share session

        self._session_id = shared_session_name
        self.close()

    def send_command(self, cmd_str, wait4reply:bool = True):
        self.connect()
        if wait4reply is False:
            self._oneway_command(cmd_str)
        else:
            self._command(cmd_str)
        self.close()


    def read_data(self, command_to_read_data:str, n_args_out: int = 1):
        self.connect()
        mat_data = self.eng.eval(str(command_to_read_data),nargout=n_args_out)

        if n_args_out > 1:
            data = []
            for val in mat_data:
                val = np.array(val)
                data.append(val.tolist())
        else:
            mat_data = np.array(mat_data)
            data = mat_data.tolist()
        self.close()

        return data
    
    def connect(self):
        self.eng = matlab.engine.connect_matlab(self._session_id)

    @Pyro4.oneway
    def close(self):
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
