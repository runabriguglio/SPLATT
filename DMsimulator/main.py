import numpy as np
import matplotlib.pyplot as plt
# from astropy.io import fits
# from scipy.sparse import csr_matrix 
# from zernike_polynomials import computeZernike as czern

from geometry import Hexagons
from finite_elements import Mesh

# from scipy.sparse.linalg import svds
from scipy.sparse.linalg import lsqr
# from scipy.sparse import csr_matrix

config_tn = '20240920'
hexes = Hexagons(config_tn)

# Show the segments mask
mask = hexes.global_mask
plt.figure()
plt.imshow(mask,origin='lower',cmap='gray')
plt.title('Segment Mask')

# Interaction matrix and initial scramble
n_modes = 11
hexes.create_interaction_matrix(n_modes)
img = hexes.segment_scramble()*10e-6
plt.figure()
plt.imshow(img,origin='lower',cmap='hot_r')
plt.title('Segment Scramble')
plt.colorbar()

w = img.data.flatten()

# # RECONSTRUCTOR
# acc = 1e-7
INT_MAT = hexes.int_mat
# command, itstop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm, var  = lsqr(INT_MAT, w, atol = acc)
# w_star = INT_MAT * command

# img_star = np.reshape(w_star, np.shape(hexes.global_mask))
# plt.figure()
# plt.imshow(np.ma.masked_array(img_star,hexes.global_mask),origin='lower',cmap='hot_r')
# plt.title('Correction command')
# plt.colorbar()
# delta_img = img - img_star
# # delta_img = np.reshape(w-w_star, np.shape(hexes.global_mask))
# masked_delta_img = np.ma.masked_array(delta_img, hexes.global_mask)

# plt.figure()
# plt.imshow(masked_delta_img,origin='lower')
# plt.title('Flattening result')
# plt.colorbar()


# # Draw hexes and inscribed circle
# hexes.draw_hex_outline()


# Mesh class
acts = Mesh(config_tn)
n_hex = len(hexes.hex_centers)
n_w = np.size(hexes.global_mask)
acts.act_influence_functions(hexes.hex_indices, hexes.local_mask, n_w, n_hex)

n_acts_per_hex = len(acts.local_act_coords)
IFF = acts.IFF

# Shape command
shape_vec = np.zeros(n_modes)
shape_vec[8] = 1
shape_cmd = np.tile(shape_vec, n_hex)
shapeW = INT_MAT * shape_cmd
# actCMD, itstop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm, var  = lsqr(IFF, shapeW, atol = 1e-5)

img = np.reshape(shapeW, np.shape(hexes.global_mask))
masked_img = np.ma.masked_array(img, hexes.global_mask)

plt.figure()
plt.imshow(masked_img, origin='lower')
plt.title('Shape command')
plt.colorbar()

# Testing IFF matrix
act_command = np.zeros(n_acts_per_hex)
command = np.tile(act_command, n_hex)
idx = np.arange(n_hex) * n_acts_per_hex + np.arange(n_hex)
command[idx] = 1
w_cmd = IFF * command
cmd_img = np.reshape(w_cmd, np.shape(hexes.global_mask))
masked_cmd = np.ma.masked_array(cmd_img, hexes.global_mask)

plt.figure()
plt.imshow(masked_cmd, origin='lower')
plt.colorbar()

# # Plotting actuator locations
# local_mesh_points = acts.local_mesh_coords
# local_act_coords = acts.local_act_coords

# plt.figure()
# plt.grid('on')
# plt.axis('equal')

# for center in hexes.hex_centers:
#     # print(center)
#     mesh_coords = local_mesh_points*hexes.hex_side_len + center
#     global_act_coords = local_act_coords*hexes.hex_side_len + center
#     # tria = Delaunay(mesh_coords) #,incremental=True
#     # plt.triplot(mesh_coords[:,0], mesh_coords[:,1], tria.simplices,c='black')
#     plt.scatter(global_act_coords[:,0], global_act_coords[:,1],c='red')
    
# plt.xlim([-2,2])
# plt.ylim([-2,2])

# Plotting local act coordinates
local_act_coords = acts.local_act_coords

plt.figure()
plt.grid('on')
plt.axis('equal')
plt.scatter(local_act_coords[:,0], local_act_coords[:,1],c='red')

# # Plotting fake influence functions
# masked_cube = acts.local_IFFs
# for img in masked_cube:
#     plt.figure()
#     plt.imshow(img,origin='lower')
    

# RECONSTRUCTOR
acc = 1e-9
command, itstop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm, var  = lsqr(IFF, w, atol = acc)
w_star = IFF * command

img_star = np.reshape(w_star, np.shape(hexes.global_mask))
plt.figure()
plt.imshow(np.ma.masked_array(img_star,hexes.global_mask),origin='lower',cmap='hot_r')
plt.title('Correction command')
plt.colorbar()
delta_img = np.reshape(w-w_star, np.shape(hexes.global_mask))
masked_delta_img = np.ma.masked_array(delta_img, hexes.global_mask)

plt.figure()
plt.imshow(masked_delta_img,origin='lower')
plt.title('Flattening result')
plt.colorbar()


