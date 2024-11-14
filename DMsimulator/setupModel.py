import numpy as np
import matplotlib.pyplot as plt

# from HexagonClass import Hexagon
from DMClass import SegmentedMirror

def dm_system_setup(TN):
    
    # Build segmented deformable mirror
    dm = SegmentedMirror(TN)
    
    # Plot global mask
    plt.figure()
    plt.imshow(dm.global_mask, origin = 'lower', cmap='gray')
    plt.title('Global Mask')
    
    # Global interaction matrix
    N_global_modes = 11
    dm.compute_global_interaction_matrix(N_global_modes)
    glob_INTMAT = dm.glob_int_mat
    tiptilt = [0,1,1,0]
    wf = glob_INTMAT * tiptilt
    dm.plot_wavefront(wf, 'Global Tip/Tilt')
    
    # Interaction matrix
    N_modes = 11
    dm.compute_interaction_matrix(N_modes)
    INTMAT = dm.int_mat
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
    
    # Global influence functions and global reconstructor
    dm.assemble_IFF_and_R_matrices()
    
    return dm



# # Actuator data
# n_acts = len(hexA.local_act_coords)
# max_x, max_y = np.shape(dm.local_mask)
# valid_len = np.sum(1-dm.local_mask)

# #Fitting error for a single segment
# k = 1
# hex_k_ids = dm.hex_valid_ids[k]
# loc_IM = INTMAT[hex_k_ids,N_modes*k:N_modes*(k+1)]
# loc_R = np.linalg.pinv(loc_IFF, rcond = 1e-15)
# rel_rms = np.zeros(N_modes)

# wf = np.zeros(np.size(hexA.local_mask))
# mask = hexA.local_mask.copy()
# fit_cube_err = []
# flat_mask = mask.flatten()
# flat_image = np.zeros(np.size(flat_mask))

# for j in np.arange(N_modes):
    
#     modes = np.zeros(N_modes)
#     modes[j] = 1 
    
#     # modal_img = loc_IM * modes
#     # w_mode = modal_img[hex_k_ids]
#     w_mode = loc_IM * modes
#     cmd = loc_R @ w_mode
#     w_cmd = loc_IFF @ cmd
#     w = w_mode - w_cmd
    
#     rel_rms[j] = np.std(w)#/np.std(w_mode)
    
#     flat_image[~flat_mask] = w
#     image = np.reshape(flat_image, np.shape(hexA.local_mask))
#     image = np.ma.masked_array(image, hexA.local_mask)
#     fit_cube_err.append(image)
#     plt.figure()
#     plt.imshow(image, origin = 'lower', cmap = 'gray')
#     plt.title('Mode ' + str(j) + ' residual\n RMS: ' + str(np.std(w)) )
#     plt.colorbar()
    
# plt.figure()
# plt.plot(rel_rms,'o')
# plt.grid('on')
# plt.title('RMS residual fitting error')

# plt.figure()
# id_err = loc_R @ loc_IFF - np.eye(len(loc_R[:,0]))
# plt.imshow(id_err, origin='lower')
# plt.title('Reconstructor error')
# plt.colorbar()


# # Plot IFF data on segments
# n_hex = len(dm.hex_valid_ids)
# glob_img = np.zeros(np.size(dm.global_mask))
# point_glob_img = np.zeros(np.size(dm.global_mask))
# for kk in range(n_hex):
#     glob_img[dm.hex_valid_ids[kk]] = hexA.sim_IFF[:,kk]#loc_IFF[:,kk]
#     # point_glob_img[dm.hex_valid_ids[kk]] = hexA.IFF[:,kk]
    
# # dm.plot_wavefront(glob_img, 'Actuator Influence Functions')
# full_img = np.reshape(glob_img, np.shape(dm.global_mask))
# img = np.ma.masked_array(full_img, dm.global_mask)
# plt.figure()
# plt.imshow(img, origin = 'lower', cmap='inferno')
# # plt.axis([1215,1620,1220,1600])
# plt.colorbar()


