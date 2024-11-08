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
n_hex = int(np.shape(INTMAT)[1]/N_modes)
cmd_ids = np.arange(N_modes-1)+1
cmd_ids = np.tile(cmd_ids,int(np.ceil(n_hex/(N_modes-1))))
cmd_ids = cmd_ids[0:n_hex]
cmd_ids = cmd_ids + N_modes*np.arange(n_hex)
modal_cmd = np.zeros(n_hex*N_modes)
modal_cmd[cmd_ids] = 1
flat_img = INTMAT * modal_cmd
dm.plot_wavefront(flat_img, 'Zernike modes')


# Testing influence functions
iff_cube = hexA.simulate_influence_functions()
IFF = hexA.IFF

n_acts = np.shape(IFF)[1]
for k in np.arange(n_acts):
    # img = iff_cube[:,:,k]
    # # img = np.ma.masked_array(img, hexA.local_mask)
    # plt.figure()
    # plt.imshow(img,origin='lower',cmap='gray')
    cmd = np.zeros(n_acts)
    cmd[k] = 1
    w = IFF*cmd
    img = np.reshape(w, np.shape(hexA.local_mask))
    img = np.ma.masked_array(img,hexA.local_mask)
    plt.figure()
    plt.imshow(img,origin='lower',cmap='gray')



# #####
# def generate_hex_mask(npix = int(2**9)):
#     L = npix/(1.+2.*COS60)
#     half_npix = int(npix/2)
    
#     # Consider points in the upper right of the hexagon
#     mask_ul = np.fromfunction(lambda i,j: i >=  np.minimum(j/COS60,L)*SIN60, [half_npix,half_npix])

#     mask = np.zeros([npix,npix], dtype=bool)

#     mask[half_npix:,0:half_npix] = mask_ul # upper right
#     mask[half_npix:,half_npix:] = np.flip(mask_ul,1) # upper left
#     mask[0:half_npix,:] = np.flip(mask[half_npix:,:],0) # lower
    
#     return mask

# hex_mask = generate_hex_mask()
# plt.figure()
# plt.imshow(hex_mask,origin='lower',cmap='gray')

# npix = int(2**9)
# coords = dm.hex_centers
# limit_coords = coords + np.max([1.+2*COS60, 2*SIN60])*dm.hex_side_len
# D = np.max(limit_coords)
# pix_coords = (coords+D)/(2*D)*npix
# mask_size = int(np.round(dm.hex_side_len*(1.+2*COS60)/(2*D)*npix))
# mask_coords = (pix_coords - mask_size*2).astype(int)

# global_mask = np.ones([npix,npix], dtype=bool)
# scaled_mask = generate_hex_mask(mask_size)
# for m_coords in mask_coords:
#     global_mask[m_coords[0]:m_coords[0]+mask_size,m_coords[1]:m_coords[1]+mask_size] = scaled_mask
  
# plt.imshow(global_mask, origin='lower', cmap = 'gray')