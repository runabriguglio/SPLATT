import numpy as np

def cw_rotate(vec, angles):
    """ Rotates vector vec of shape [2,n_points]
    clockwise about [0,0] by all angles in input
    Outputs a vector of shape vec [n_points,2]"""
    
    L = len(angles)
        
    n_pts = int(np.size(vec)/2)
    rot_vec = np.zeros([2,n_pts*L])
    
    for k, angle in enumerate(angles):
        c = np.cos(angle)
        s = np.sin(angle)
        rot_mat = [[c,s],[-s,c]]
        
        if np.shape(vec)[0] != 2:
            aux_vec = rot_mat @ vec.transpose()
        else:
            aux_vec = rot_mat @ vec
            
        rot_vec[:,k*n_pts:(k+1)*n_pts] = np.reshape(aux_vec,[2,n_pts])
    
    return rot_vec.T


def rotate_by_60deg(vec):
    # Wrapper to cw_rotate()
    angles = np.array([np.pi/3.])
    rot_vec = cw_rotate(vec, angles)
    
    return rot_vec
