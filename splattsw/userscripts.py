from m4.mini_OTT import timehistory as th
import numpy as np
from splattsw.splatt_utilities import read_buffer_data

def analyze_buf_step(tn, modes, bits2m:float = 2**-26, bits2N:float = 1.0/99643):
    data, time_vec = read_buffer_data(tn)

    pos = data['sabi32_Distance']*bits2m
    cur = data['sabi32_pidCoilOut']*bits2N
    pos_cmd = data['sabu16_position']*bits2m
    cur_cmd = data['sabi16_force']*bits2N

    pos_modes = get_modal_projection(pos, modes)
    cur_modes = get_modal_projection(cur, modes)
    pos_cmd_modes = get_modal_projection(pos_cmd, modes)
    cur_cmd_modes = get_modal_projection(cur_cmd, modes)

    no_rigid_modes = pos_modes[4:,:]
    modes_rms = no_rigid_modes.std(axis=0)

    data_modes = {
        'position' : pos_modes.T,
        'force' : cur_modes.T,
        'position_command' : pos_cmd_modes.T,
        'force_command' : cur_cmd_modes.T,
        'pos_rms' : modes_rms.T,
        'time' : time_vec.T}

    return data_modes


def get_modal_projection(data, modes):
    delta_data = data - np.reshape(data[:,0],[19,1])
    data_modes = modes.T @ delta_data

    return data_modes


def analyze_opt_step(tn):
    fl = th.fileList(tn)
    cc = th.cubeFromList(fl)
    dd = cc-cc[0]
    zz = []
    for i in dd:
        zz.append(th.removeZernike(i,[1,2,3,4]))
    zz = np.ma.masked_array(zz)
    vv = zz.std( axis=(1,2))
    return vv

