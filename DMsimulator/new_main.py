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
full_img = np.reshape(flat_img, np.shape(dm.global_mask))
img = np.ma.masked_array(full_img, dm.global_mask)

plt.figure()
plt.imshow(img, origin = 'lower', cmap='winter')
plt.title('Global Power')

glob_tilt = [0,1,-1,0]
flat_img = global_INTMAT * glob_tilt
full_img = np.reshape(flat_img, np.shape(dm.global_mask))
img = np.ma.masked_array(full_img, dm.global_mask)

plt.figure()
plt.imshow(img, origin = 'lower', cmap='winter')
plt.title('Global Tip/Tilt')

# Initial segment scramble
scrambled_img = dm.segment_scramble()
plt.figure()
plt.imshow(scrambled_img, origin = 'lower', cmap = 'hot')
plt.title('Segment scramble')

# Testing single segment commands
N = 8
local_INTMAT_N = INTMAT[:,N_modes*(N-1):N_modes*N]
cmd = np.zeros(N_modes)
cmd[8] = 1
flat_img = local_INTMAT_N * cmd
full_img = np.reshape(flat_img, np.shape(dm.global_mask))
img = np.ma.masked_array(full_img, dm.global_mask)

plt.figure()
plt.imshow(img, origin = 'lower', cmap='winter')
plt.title('Local Trefoil')