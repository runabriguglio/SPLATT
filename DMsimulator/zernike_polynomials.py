import numpy as np

def computeZernike(noll_number, mask):

    X,Y = np.shape(mask)
    
    match noll_number:
        case 1: # Piston
            mode = np.ones([X,Y])
        case 2: # Tip
            mode = np.fromfunction(lambda i,j: 2.*j/Y - 1., [X,Y])
        case 3: # Tilt
            mode = np.fromfunction(lambda i,j: 2.*i/X - 1., [X,Y])
        case 4: # Focus
            mode = np.fromfunction(lambda i,j: ((j-Y/2.)/Y)**2+((i-X/2.)/Y)**2, [X,Y])
        case _:
            print("Zernike mode too high: not yet implemented")
            return

    masked_mode = np.ma.masked_array(mode,mask)

    masked_data = masked_mode.data[~masked_mode.mask]

    # Normalization: null mean and unit STD
    if noll_number > 1:
        masked_data = (masked_data - np.mean(masked_data))/np.std(masked_data)
        
    return masked_data