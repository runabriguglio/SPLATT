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

    def send(self, cmd_str:str):
        self.connect()
        
        try:
            self.eng.eval(str(cmd_str)+';',nargout=0)
        except:
            self.close()
            raise Pyro4.errors.DaemonError(f"Execution of matlab command '{cmd_str}' encountered an error")

        self.close()


    def read(self, cmd_returning_data:str, n_args_out: int = 1):
        self.connect()

        try:
            mat_data = self.eng.eval(str(cmd_returning_data),nargout=n_args_out)
        except:
            self.close()
            raise Pyro4.errors.DaemonError(f"Execution of matlab command '{cmd_returning_data}' encountered an error")

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

    @Pyro4.oneway
    def oneway_send(self, cmd_str:str):
        self.send(cmd_str)
    
    def connect(self):
        self.eng = matlab.engine.connect_matlab(self._session_id)

    @Pyro4.oneway
    def close(self):
        self.eng.quit()


def main():

    matlab_eng = MatlabEngine()

    Pyro4.Daemon.serveSimple( {matlab_eng: "matlab_engine"},
             host="193.206.155.220", port=9090, ns=False, verbose=True)


if __name__=="__main__":
    main()


    
