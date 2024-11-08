import numpy as np
import matplotlib.pyplot as plt
# from scipy.interpolate import griddata
# from scipy.sparse.linalg import lsqr

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


# Initial segment scramble
scrambled_img = dm.segment_scramble()
w_scramble = scrambled_img.flatten()
dm.plot_wavefront(w_scramble, 'Segment scramble')

# Testing influence functions
hexA.simulate_influence_functions()
loc_IFF = hexA.IFF

# pinv_IFF = np.linalg.pinv(local_IFF)
# flattone = w_scramble.copy()

# for k in np.arange(n_hex):
#     hex_k_ids = dm.hex_valid_ids[k]
#     w = w_scramble[hex_k_ids]
#     acc = 1e-9
#     #command, itstop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm, var  = lsqr(local_IFF, w, atol = acc)
#     command = pinv_IFF @ w
#     w_star = local_IFF @ command
    
#     flattone[hex_k_ids] -= w_star
#     dm.plot_wavefront(flattone, 'Hex ' + str(k+1) + ' correction')

# flat_img = np.reshape(flattone,np.shape(dm.global_mask))
# flat_img = np.ma.masked_array(flat_img, dm.global_mask)
# plt.figure()
# fig=plt.imshow(flat_img,origin='lower',cmap='gray')
# cmap=plt.colorbar()
# masked_img_data = flat_img.data[~flat_img.mask]
# fig.set_clim([min(masked_img_data),max(masked_img_data)])
# plt.axis([950,1000,950,1000])
# flat_rel_rms = np.std(masked_img_data)/np.std(w_scramble)

#Fitting error for a single segment
k = 1
hex_k_ids = dm.hex_valid_ids[k]
loc_IM = INTMAT[hex_k_ids,N_modes*k:N_modes*(k+1)]
loc_R = np.linalg.pinv(loc_IFF, rcond = 1e-15)
rel_rms = np.zeros(N_modes)

wf = np.zeros(np.size(dm.global_mask))

for j in np.arange(N_modes-1):
    
    modes = np.zeros(N_modes)
    modes[j+1] = 1 
    
    # modal_img = loc_IM * modes
    # w_mode = modal_img[hex_k_ids]
    w_mode = loc_IM * modes
    cmd = loc_R @ w_mode
    w_cmd = loc_IFF @ cmd
    w = w_mode - w_cmd
    
    rel_rms[j+1] = np.std(w)/np.std(w_mode)
    
    wf[hex_k_ids] = w
    img = np.reshape(wf, np.shape(dm.global_mask))
    img = np.ma.masked_array(img, dm.global_mask)
    plt.figure()
    plt.imshow(img, cmap='gray', origin = 'lower')
    plt.title('Mode ' + str(j+1) + ' residual\n RMS: ' + str(np.std(w)) )
    plt.colorbar()
    plt.axis([1220,1630,1225,1600])
    

def meters_to_pixels(coords, mask = hexA.local_mask, pixel_scale = hexA.pix_scale):
    
    x, y = coords[:,0], coords[:,1]
    max_x, max_y = np.shape(mask)
    
    pix_coords = np.zeros_like(coords)
    pix_coords[:,0] = (x * pixel_scale + max_y/2).astype(int)
    pix_coords[:,1] = (y * pixel_scale + max_x/2).astype(int)
    
    return pix_coords
    
pix_act_coords = meters_to_pixels(hexA.local_act_coords)
plt.figure()
plt.scatter(pix_act_coords[:,0],pix_act_coords[:,1],s=100)
plt.axis('equal')

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
