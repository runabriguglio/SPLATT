import numpy as np

def computeZernike(noll_number, rect_shape):
    Y = rect_shape[1]
    X = rect_shape[0]
    
    match noll_number:
        case 0: # Piston
            norm_mode = np.ones([X,Y])
        case 1: # Tip
            mode = np.fromfunction(lambda i,j: 2.*j/Y - 1., [X,Y])
        case 2: # Tilt
            mode = np.fromfunction(lambda i,j: 2.*i/X - 1., [X,Y])
        case 3: # Focus
            mode = np.fromfunction(lambda i,j: ((j-Y/2.)/Y)**2+((i-X/2.)/Y)**2, [X,Y])
        case _:
            print("Zernike mode too high: not yet implemented")
            return
        
    # Normalization: null mean and unit STD
    if noll_number > 0:
        norm_mode = (mode - np.mean(mode))/np.std(mode)
        
    return norm_mode