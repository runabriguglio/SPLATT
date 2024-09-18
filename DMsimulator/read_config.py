import configparser
import numpy as np

config = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))

def readConfig(path):
    config.read(path)

    dm_par = config['DM']
    opt_par = config['OPT']
    
    g = dm_par['hex_gap']
    l_hex = dm_par['hex_side']
    n_rings = dm_par['n_rings']

    ang = opt_par['cw_rot_angle']
    pix_scale = opt_par['pixel_scale']

    par = np.array([g,l_hex,n_rings,pix_scale,ang])
    par = par.astype(float)

    return par
