import numpy as np

def cw_rotate(vec, angle):
    c = np.cos(angle)
    s = np.sin(angle)
    rot_mat = [[c,s],[-s,c]]
    
    if np.shape(vec)[0] > 2:
        aux_vec = rot_mat @ vec.transpose()
        rot_vec = aux_vec.transpose()
    else:
        rot_vec = rot_mat @ vec
    
    return rot_vec


def rotate_by_60deg(vec):
    # Wrapper to cw_rotate()
    cw_angle = np.pi/3.
    rot_vec = cw_rotate(vec, cw_angle)
    
    return rot_vec
