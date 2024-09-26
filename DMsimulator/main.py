import numpy as np
import matplotlib.pyplot as plt
from geometry import Hexagons
# from astropy.io import fits
# from scipy.sparse import csr_matrix 
# from zernike_polynomials import computeZernike as czern

# Mesh testing
from finite_elements import Mesh
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



# Mesh class
acts = Mesh(config_tn)
acts.act_influence_functions(hexes.local_mask)


# Plotting actuator locations
local_mesh_points = acts.local_mesh_coords
local_act_coords = acts.local_act_coords

plt.figure()
plt.grid('on')
plt.axis('equal')

for center in hexes.hex_centers:
    # print(center)
    mesh_coords = local_mesh_points*hexes.hex_side_len + center
    global_act_coords = local_act_coords*hexes.hex_side_len + center
    tria = Delaunay(mesh_coords) #,incremental=True
    plt.triplot(mesh_coords[:,0], mesh_coords[:,1], tria.simplices,c='black')
    plt.scatter(global_act_coords[:,0], global_act_coords[:,1],c='red')
    
plt.xlim([-2,2])
plt.ylim([-2,2])


# Plotting fake influence functions
masked_cube = acts.local_IFFs
for img in masked_cube:
    plt.figure()
    plt.imshow(img,origin='lower')



