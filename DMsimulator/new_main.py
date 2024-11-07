import numpy as np
import matplotlib.pyplot as plt
# from scipy.interpolate import griddata

from HexagonClass import Hexagon
from DMClass import DM

config_tn = '20240920'

hexA = Hexagon(config_tn)
dm = DM(config_tn)

plt.figure()
plt.imshow(dm.global_mask, origin = 'lower', cmap='gray')
plt.title('Global Mask')

# Interaction matrices
N_modes = 11
dm.compute_interaction_matrix(N_modes)
INTMAT = dm.int_mat

N_global_modes = 4
dm.compute_global_interaction_matrix(N_global_modes)
global_INTMAT = dm.glob_int_mat

glob_pow = [0,0,0,-1]
flat_img = global_INTMAT * glob_pow
dm.plot_wavefront(flat_img, 'Global Power')

glob_tilt = [0,1,-1,0]
flat_img = global_INTMAT * glob_tilt
dm.plot_wavefront(flat_img, 'Global Tip/Tilt')

# Initial segment scramble
scrambled_img = dm.segment_scramble()
plt.figure()
plt.imshow(scrambled_img, origin = 'lower', cmap = 'inferno_r')
plt.title('Segment scramble')

# Testing single segment commands
n_hexes = int(np.shape(INTMAT)[1]/N_modes)
cmd_ids = np.arange(N_modes-1)+1
cmd_ids = np.tile(cmd_ids,int(np.ceil(n_hexes/(N_modes-1))))
cmd_ids = cmd_ids[0:n_hexes]
cmd_ids = cmd_ids + N_modes*np.arange(n_hexes)
modal_cmd = np.zeros(n_hexes*N_modes)
modal_cmd[cmd_ids] = 1
flat_img = INTMAT * modal_cmd
dm.plot_wavefront(flat_img, 'Zernike modes')


# Testing influence functions
hexA.simulate_influence_functions()
IFF = hexA.IFF