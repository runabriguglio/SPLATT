import numpy as np
import matplotlib.pyplot as plt
from geometry import Hexagons

# from scipy.sparse import csr_matrix 
# from zernike_polynomials import computeZernike as czern

configfile = 'config_parameters.yaml'
hexes = Hexagons(configfile)

hexes.find_hex_coordinates()
# plt.plot(c_coords[:,0],c_coords[:,1],'o')

hexes.define_hex_mask()
# mask = hexes.hex_mask
# plt.imshow(mask,origin = 'lower')

hexes.define_segmented_mask()
mask = hexes.full_mask
plt.imshow(mask,origin='lower')

n_modes = 3
hexes.calculate_interaction_matrix(n_modes)
img = hexes.segment_scramble()
plt.figure()
plt.imshow(img,origin='lower')


# tip = czern(2,hexes.hex_mask)
# tip1 = np.zeros(np.size(hexes.full_mask))
# tip1[hexes.hex_indices[0,:]] = tip
# tip1 = np.reshape(tip1,hexes.full_mask.shape)
# plt.figure()
# plt.imshow(tip1,origin='lower')

# tilt = czern(3,hexes.hex_mask)
# tilt8 = np.zeros(np.size(hexes.full_mask))
# tilt8[hexes.hex_indices[7,:]] = tilt
# tilt8 = np.reshape(tilt8,hexes.full_mask.shape)
# plt.figure()
# plt.imshow(tilt8,origin='lower')

# mat1 = hexes.int_mat[:,6:18]
# flat_img = mat1 @ np.tile(np.array([0.,0.,1.]),4) #mat1.dot(np.array([0.,1.,0.]))
# img1 = np.reshape(flat_img, hexes.full_mask.shape)
# plt.figure()
# plt.imshow(img1,origin='lower')

# flat_indices = hexes.hex_indices.flatten()
# mode_indices = np.arange(n_modes)
# X,Y = np.meshgrid(flat_indices,mode_indices)



# hex_mask_len = hexes.hex_mask.size
# flat_mask = hexes.hex_mask.flatten()

# [Lx,Ly] = hexes.hex_mask.shape
# tip = np.fromfunction(lambda i,j: 2.*j - 1., [Lx,Ly])
# tip = (tip-np.mean(tip))/np.std(tip)
# flat_tip = tip.flatten()

# # Masked array
# masked_ar = np.ma.array(flat_tip,mask=flat_mask)
# masked_img = np.reshape(masked_ar,hexes.hex_mask.shape)
# plt.figure()
# plt.imshow(masked_img,origin='lower')

# # Flat array, using indices
# indices = np.arange(hex_mask_len)
# valid_indices = indices[~flat_mask]
# masked_flat_tip = np.zeros_like(flat_tip)
# masked_flat_tip[valid_indices] = flat_tip[valid_indices]
# masked_img = np.reshape(masked_flat_tip,hexes.hex_mask.shape)
# masked_img = np.ma.masked_array(masked_img, hexes.hex_mask)
# plt.figure()
# plt.imshow(masked_img,origin='lower')

# # Sparse matrix
# x = np.arange(Lx)
# y = np.arange(Ly)
# X,Y = np.meshgrid(y,x)
# valid_X = X[~hexes.hex_mask]
# valid_Y = Y[~hexes.hex_mask]

# data = flat_tip[valid_indices]
# row = valid_Y
# col = valid_X
# sparseMatrix = csr_matrix((data, (row, col)),  
#                           shape = hexes.hex_mask.shape).toarray() 
# # plt.figure()
# # plt.imshow(sparseMatrix,origin='lower')
# masked_img = np.ma.masked_array(sparseMatrix, hexes.hex_mask)
# plt.figure()
# plt.imshow(masked_img,origin='lower')

