import numpy as np
# import matplotlib.pyplot as plt

def compute_zernike(noll_number, mask):
    """ Defines the Zernike polynomials identified by the Noll number in input
    on a mask (boolean array). The polynomials are computed in a circle
    inscribed in the mask, and the masked data is then normalized on the 
    valid area as follows:
    data = ma.data[~ma.mask], data = (data - mean(data))/std(data)"""

    X,Y = np.shape(mask)
    
    # Determine circle radius on to which define the Zernike
    r = np.max([X,Y])/2
    theta = lambda i,j: np.arctan2((j-Y/2.)/r,(i-X/2.)/r)
    rho = lambda i,j: np.sqrt(((j-Y/2.)/r)**2+((i-X/2.)/r)**2)
    
    match noll_number:
        case 1: # Piston
            mode = np.ones([X,Y])
        case 2: # Tip
            mode = np.fromfunction(lambda i,j: 2.*j/Y - 1., [X,Y])
        case 3: # Tilt
            mode = np.fromfunction(lambda i,j: 2.*i/X - 1., [X,Y])
        case 4: # Focus
            mode = np.fromfunction(lambda i,j: ((j-Y/2.)/Y)**2+((i-X/2.)/Y)**2, [X,Y])
        case 5: # Vertical astigmatism
            mode = np.fromfunction(lambda i,j: (rho(i,j)**2)*np.cos(2.*theta(i,j)), [X,Y]) #np.cos(2*np.atan2((j-Y/2.)/(i-X/2.)))*((j-Y/2.))**2+((i-X/2.))**2
        case 6: # Oblique astigmatism
            mode = np.fromfunction(lambda i,j: (rho(i,j)**2)*np.sin(2.*theta(i,j)), [X,Y])
        case 7: # Vertical coma
            mode = np.fromfunction(lambda i,j: (3.*rho(i,j)**3-2.*rho(i,j))*np.sin(theta(i,j)), [X,Y])
        case 8: # Horizontal coma
            mode = np.fromfunction(lambda i,j: (3.*rho(i,j)**3-2.*rho(i,j))*np.cos(theta(i,j)), [X,Y])
        case 9: # Vertical trefoil
            mode = np.fromfunction(lambda i,j: (rho(i,j)**3)*np.sin(3.*theta(i,j)), [X,Y])
        case 10: # Oblique trefoil
            mode = np.fromfunction(lambda i,j: (rho(i,j)**3)*np.cos(3.*theta(i,j)), [X,Y])
        case 11: # Primary spherical
            mode = np.fromfunction(lambda i,j: (rho(i,j)**4-rho(i,j)**2), [X,Y])
        case _:
            print("Zernike mode too high: not yet implemented")
            return

    masked_mode = np.ma.masked_array(mode,mask)

    masked_data = masked_mode.data[~masked_mode.mask]
    
    # masked_mode = (masked_mode - np.mean(masked_data))/np.std(masked_data)
    # plt.figure()
    # plt.imshow(masked_mode,origin='lower')

    # Normalization: null mean and unit STD
    if noll_number > 1:
        masked_data = (masked_data - np.mean(masked_data))/np.std(masked_data)
        
    return masked_data