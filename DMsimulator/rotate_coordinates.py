import numpy as np

def cw_rotate(vec, angles):
    
    n_pts = int(np.size(vec)/2)
    res_vec = np.zeros([n_pts*len(angles),2])
    
    for k, angle in enumerate(angles):
        c = np.cos(angle)
        s = np.sin(angle)
        rot_mat = [[c,s],[-s,c]]
        
        if n_pts > 2:
            aux_vec = rot_mat @ vec.transpose()
            rot_vec = aux_vec.transpose()
        else:
            rot_vec = rot_mat @ vec
            
        res_vec[k*n_pts:(k+1)*n_pts,:] = rot_vec
    
    return res_vec


def rotate_by_60deg(vec):
    # Wrapper to cw_rotate()
    angles = np.array([np.pi/3.])
    rot_vec = cw_rotate(vec, angles)
    
    return rot_vec
