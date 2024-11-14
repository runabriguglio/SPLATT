import configparser
import numpy as np

config = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))

def read_config(path):
    config.read(path)

     # DM configuration parameters
    dm_conf = config['DM']
    g = dm_conf['hex_gap']
    l_hex = dm_conf['hex_side']
    n_rings = dm_conf['n_rings']
    act_pitch = dm_conf['act_pitch']
    act_r = dm_conf['act_radius']

    dm_par = np.array([g, l_hex, n_rings, act_pitch, act_r])
    dm_par = dm_par.astype(float)

    # Optical configuration parameters
    opt_conf = config['OPT']
    ang = opt_conf['cw_rot_angle']
    pix_scale = opt_conf['pixel_scale']
    pup_x = opt_conf['pupil_x']
    pup_y = opt_conf['pupil_y']
    opt_rad = opt_conf['opt_radius']

    opt_par = np.array([pix_scale, ang, pup_x, pup_y, opt_rad])
    opt_par = opt_par.astype(float)

    return dm_par, opt_par
