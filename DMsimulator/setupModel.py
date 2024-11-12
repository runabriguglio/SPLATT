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

# Actuator data
n_acts = len(hexA.local_act_coords)
max_x, max_y = np.shape(hexA.local_mask)
valid_len = np.sum(1-hexA.local_mask)

#Fitting error for a single segment
k = 1
hex_k_ids = dm.hex_valid_ids[k]
loc_IM = INTMAT[hex_k_ids,N_modes*k:N_modes*(k+1)]
loc_R = np.linalg.pinv(loc_IFF, rcond = 1e-15)
rel_rms = np.zeros(N_modes)

wf = np.zeros(np.size(dm.global_mask))

for j in np.arange(N_modes):
    
    modes = np.zeros(N_modes)
    modes[j] = 1 
    
    # modal_img = loc_IM * modes
    # w_mode = modal_img[hex_k_ids]
    w_mode = loc_IM * modes
    cmd = loc_R @ w_mode
    cmd[n_acts:] = 0*cmd[n_acts:]
    w_cmd = loc_IFF @ cmd
    w = w_mode - w_cmd
    
    rel_rms[j] = np.std(w)#/np.std(w_mode)
    
    wf[hex_k_ids] = w
    img = np.reshape(wf, np.shape(dm.global_mask))
    img = np.ma.masked_array(img, dm.global_mask)
    plt.figure()
    plt.imshow(img, cmap='gray', origin = 'lower')
    plt.title('Mode ' + str(j) + ' residual\n RMS: ' + str(np.std(w)) )
    plt.colorbar()
    plt.axis([1200,1600,1200,1600])
    
plt.figure()
plt.plot(rel_rms,'o')
plt.grid('on')
plt.title('RMS residual fitting error')

plt.figure()
id_err = loc_R @ loc_IFF - np.eye(len(loc_R[:,0]))
plt.imshow(id_err, origin='lower')
plt.title('Reconstructor error')
plt.colorbar()

def meters_to_pixels(coords, mask = hexA.local_mask, pixel_scale = hexA.pix_scale):
    
    x, y = coords[:,0], coords[:,1]
    max_x, max_y = np.shape(mask)
    
    pix_coords = np.zeros_like(coords)
    pix_coords[:,0] = (x * pixel_scale + max_y/2).astype(int)
    pix_coords[:,1] = (y * pixel_scale + max_x/2).astype(int)
    
    return pix_coords

from rotate_coordinates import rotate_by_60deg as rot60
from rotate_coordinates import cw_rotate

# Useful variables
SIN60 = np.sin(np.pi/3.)
COS60 = np.cos(np.pi/3.)

# 'Ghost' actuators
k = int((np.sqrt(12*n_acts-3)-3)/6)
pit = hexA.act_pitch/hexA.hex_len
rad = hexA.act_radius/hexA.hex_len
dx = 2*rad+pit
y = np.linspace(SIN60*k*dx,0.,k+1)
x = (k+1)*dx - COS60/SIN60 * y
n = k+1
p = np.zeros([6*n,2])
p[0:n,0] = x
p[0:n,1] = y
p[n:2*n,:] = rot60(p[0:n,:].T)
p[2*n:3*n,:] = rot60(p[n:2*n,:].T)
p[3*n:,:] = cw_rotate(p[0:3*n,:].T,np.array([np.pi]))
ghost_act_coords = p*hexA.hex_len
n_ghost_acts = len(ghost_act_coords)
    
pix_act_coords = meters_to_pixels(hexA.local_act_coords)
ghost_pix_act_coords = meters_to_pixels(ghost_act_coords)
plt.figure()
plt.scatter(pix_act_coords[:,0],pix_act_coords[:,1],s=100)
plt.scatter(ghost_pix_act_coords[:,0],ghost_pix_act_coords[:,1],s=100)
plt.axis([0,400,0,346])

# R = hexA.act_radius*hexA.pix_scale

# act_mask = np.zeros_like(hexA.local_mask)
# iff_cube = np.zeros([max_x,max_y, len(pix_act_coords)])
    
# for k in range(len(pix_act_coords)):
#     act_coord = pix_act_coords[k,:]
#     act_iff = np.fromfunction(lambda i,j: (i-act_coord[1])**2 + (j-act_coord[0])**2 < R**2, [max_x,max_y])
#     # plt.figure()
#     act_mask += act_iff
#     # img = np.ma.masked_array(act_iff,hexA.local_mask)
#     # plt.imshow(act_iff,origin='lower')
#     iff_cube[:,:,k]= np.ones_like(hexA.local_mask)*act_iff
    
# # plt.figure()
# # plt.imshow(act_mask,origin='lower',cmap='gray')    

# cube_mask = np.tile(act_mask,len(pix_act_coords))
# cube_mask= np.reshape(cube_mask,np.shape(iff_cube),order='F')
# iff_cube = np.ma.masked_array(iff_cube,~cube_mask)

# # for k in range(len(pix_act_coords)):
# #     plt.figure()
# #     plt.imshow(iff_cube[:,:,k],origin='lower' ) 
    
# act_x_pix = np.arange(max_x)
# act_y_pix = np.arange(max_y)

# X,Y = np.meshgrid(act_y_pix,act_x_pix)
# xx = X[act_mask] #valid_X
# yy = Y[act_mask] #valid_Y

# # npix_x, npix_y = np.shape(hexA.local_mask)  
# new_x = np.linspace(min(xx), max(xx), max_y) # x coordinate is the img column!
# new_y = np.linspace(min(yy), max(yy), max_x) # y coordinate is the img row!
# gx, gy = np.meshgrid(new_x, new_y)

# from scipy.interpolate import griddata
# act_data = iff_cube.data[~iff_cube.mask]
# act_data = np.reshape(act_data, [len(xx),n_acts])
# img_cube = griddata((xx, yy), act_data, (gx, gy), fill_value = 0., method = 'cubic')

pix_coords = np.zeros([max_x*max_y,2])
pix_coords[:,0] = np.repeat(np.arange(max_x),max_y)
pix_coords[:,1] = np.tile(np.arange(max_y),max_x)

from scipy.interpolate import CloughTocher2DInterpolator as interp2D

real_pix_coords = np.zeros([n_acts+n_ghost_acts,2])
real_pix_coords[:n_acts,0] = pix_act_coords[:,1]
real_pix_coords[:n_acts,1] = pix_act_coords[:,0]
real_pix_coords[n_acts:,0] = ghost_pix_act_coords[:,1]
real_pix_coords[n_acts:,1] = ghost_pix_act_coords[:,0]

N = n_acts+n_ghost_acts
img_cube = np.zeros([max_x,max_y,N])

for k in range(N):
    act_data = np.zeros(n_acts+n_ghost_acts)
    act_data[k] = 1
    interpolator = interp2D(real_pix_coords, act_data, tol = 1e-8)
    flat_img = interpolator(pix_coords)
    img_cube[:,:,k] = np.reshape(flat_img, [max_x,max_y])


# Masked array
cube_mask = np.tile(hexA.local_mask,N)
cube_mask = np.reshape(cube_mask, np.shape(img_cube), order = 'F')
cube = np.ma.masked_array(img_cube,cube_mask)

# Save valid data to IFF (full) matrix
flat_cube = cube.data[~cube.mask]
local_IFF = np.reshape(flat_cube, [valid_len, N])

loc_IFF = np.array(local_IFF)

# Plot IFF data on segments
n_hex = len(dm.hex_valid_ids)
glob_img = np.zeros(np.size(dm.global_mask))
point_glob_img = np.zeros(np.size(dm.global_mask))
for kk in range(n_hex):
    glob_img[dm.hex_valid_ids[kk]] = loc_IFF[:,kk]
    point_glob_img[dm.hex_valid_ids[kk]] = hexA.IFF[:,kk]
    
# dm.plot_wavefront(glob_img, 'Actuator Influence Functions')
full_img = np.reshape(glob_img, np.shape(dm.global_mask))
img = np.ma.masked_array(full_img, dm.global_mask)
plt.figure()
plt.imshow(img, origin = 'lower', cmap='inferno')
# plt.axis([1215,1620,1220,1600])
plt.colorbar()

# # dm.plot_wavefront(point_glob_img, 'Point-like Influence Functions')
# full_img = np.reshape(point_glob_img, np.shape(dm.global_mask))
# img = np.ma.masked_array(full_img, dm.global_mask)
# plt.figure()
# plt.imshow(img, origin = 'lower', cmap='inferno')
# plt.axis([1200,1600,1200,1600])
# plt.colorbar()

# # Project wavefront on mask
# mask = hexA.local_mask.copy()
# flat_mask = mask.flatten()
# image = np.zeros(np.size(flat_mask))
# image[~flat_mask] = loc_IFF[:,15].flatten()
# image = np.reshape(image, np.shape(hexA.local_mask))
# image = np.ma.masked_array(image, hexA.local_mask)

# # Plot image
# plt.figure()
# plt.imshow(image, origin = 'lower', cmap = 'inferno')
# plt.colorbar()

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
