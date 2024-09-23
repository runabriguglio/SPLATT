import numpy as np
import matplotlib.pyplot as plt
from geometry import Hexagons
# from astropy.io import fits
# from scipy.sparse import csr_matrix 
# from zernike_polynomials import computeZernike as czern

from scipy.spatial import Delaunay

config_tn = '20240920'
hexes = Hexagons(config_tn)

mask = hexes.full_mask
plt.figure()
plt.imshow(mask,origin='lower',cmap='gray')

n_modes = 11
hexes.calculate_interaction_matrix(n_modes)
img = hexes.segment_scramble()
plt.figure()
# plt.grid('on')
plt.imshow(img,origin='lower',cmap='hot')

hexes.draw_hex_outline()

''''''''
SIN60 = np.sin(np.pi/3.)
COS60 = np.cos(np.pi/3.)
hex_sides = np.zeros([7,2])
hex_sides[0,:] = np.array([-2*COS60, 0.])
hex_sides[1,:] = np.array([-0.5, SIN60])
hex_sides[2,:] = np.array([-hex_sides[1,0],hex_sides[1,1]])
hex_sides[3,:] = -hex_sides[0,:]
hex_sides[4,:] = -hex_sides[1,:]
hex_sides[5,:] = -hex_sides[2,:]
hex_sides[6,:] = hex_sides[0,:]

len_vec = np.arange(0.1,1.2,0.05)
rep_len_vec = np.repeat(len_vec,hex_sides.size)
rep_hex_sides = np.tile(hex_sides.flatten(),len(len_vec))
points = rep_hex_sides*rep_len_vec
points = np.reshape(points,[2,int(len(points)/2)])

tria = Delaunay(hex_sides)#,incremental=True)
plt.figure()
plt.triplot(hex_sides[:,0], hex_sides[:,1], tria.simplices)
plt.plot(hex_sides[:,0], hex_sides[:,1], 'o')
plt.show()

# astig1 = czern(5,hexes.hex_mask)
# astig2 = czern(6,hexes.hex_mask)
# coma1 = czern(7,hexes.hex_mask)
# coma2 = czern(8,hexes.hex_mask)
# tref1 = czern(9,hexes.hex_mask)
# tref2 = czern(10,hexes.hex_mask)
# spher = czern(11,hexes.hex_mask)




