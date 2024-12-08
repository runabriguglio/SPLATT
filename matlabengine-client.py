import Pyro4

eng = Pyro4.Proxy("PYRO:matlab_engine@193.206.155.220:9090")
print('Connected! Starting the MATLAB engine...')
eng.start_engine()
print("MATLAB engine started! Remember to quit when you are done with 'eng.stop_engine()'")


# class MatlabEngineClient():
    
#     def __init__(self):
#         self.eng = Pyro4.Proxy("PYRO:matlab_engine@193.206.155.220:9090")
#         print('Connected! Starting the MATLAB engine...')
#         self.eng.start_engine()
#         print("MATLAB engine started! Remember to quit when you are done with 'eng.stop_engine()'")
