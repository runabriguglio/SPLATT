import numpy as np
import matplotlib.pyplot as plt
from geometry import Hexagons
# from astropy.io import fits
# from scipy.sparse import csr_matrix 
# from zernike_polynomials import computeZernike as czern

# Mesh testing
from finite_elements import segment_act_coordinates
from finite_elements import mesh_from_act_coordinates
from scipy.spatial import Delaunay

config_tn = '20240920'
hexes = Hexagons(config_tn)

# Show the segments mask
mask = hexes.global_mask
plt.figure()
plt.imshow(mask,origin='lower',cmap='gray')

# Interaction matrix and initial scramble
n_modes = 11
hexes.calculate_interaction_matrix(n_modes)
img = hexes.segment_scramble()
plt.figure()
# plt.grid('on')
plt.imshow(img,origin='lower',cmap='hot')

# Draw hexes and inscribed circle
hexes.draw_hex_outline()


# Testing mesh
hex_side_len = hexes.hex_side_len
act_pitch = 0.3
capsens_radius = 0.04
points_per_side = 9 # very low for visualization purposes 

local_act_coords = segment_act_coordinates(act_pitch, capsens_radius, hex_side_len)
points = mesh_from_act_coordinates(points_per_side,local_act_coords, capsens_radius, hex_side_len)

act_size = capsens_radius * 250 * np.ones(len(local_act_coords)) # plotting only

plt.figure()
plt.grid('on')
plt.axis('equal')

for center in hexes.hex_centers:
    # print(center)
    coords = points*hex_side_len + center
    global_act_coords = local_act_coords*hex_side_len + center
    tria = Delaunay(coords) #,incremental=True
    plt.triplot(coords[:,0], coords[:,1], tria.simplices,c='black')
    plt.scatter(global_act_coords[:,0], global_act_coords[:,1],s=act_size,c='red')

plt.xlim([-2,2])
plt.ylim([-2,2])




