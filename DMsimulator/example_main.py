import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

import utilities as utils

config_tn = '20240920'
# config_tn = 'JWST'
dm = utils.dm_system_setup(config_tn)

# Segment Scramble
utils.segment_scramble(dm, apply_shape=True)

# # Fitting errors
# seg0 = dm.segment[0]
# local_fit_err = utils.fitting_error_plots(seg0.mask, seg0.ZM, seg0.IFF, seg0.R)

# global_fit_err = utils.fitting_error_plots(dm.mask, dm.glob_ZM, dm.IFF, dm.R, dm.valid_ids.flatten())

# Test triangulation
#from scipy.spatial import Delaunay
import matrix_calculator as matcalc
#import meshio as mymesh

coords = dm.segment[0].act_coords
pps = 30
act_radius = 0.002

points = matcalc._define_mesh(coords, pps, act_radius)

np.savetxt('../points.txt',points,fmt='%1.6f',delimiter=',')

#tria = Delaunay(points, furthest_site=False)
#ids = np.arange(len(tria.simplices))
#for k in range(len(tria.simplices)):
#    p1,p2 = tria.points[tria.simplices[k,1:3]]
#    dx = p2[0] - p1[0]
#    dy = p2[1] - p1[1]
#    if np.sqrt(dx**2+dy**2) > 0.2:
#        ids[k] = 0
#ids = ids[ids > 0]
#tria.simplices=tria.simplices[ids,:]
#cells = [("triangle", tria.simplices[:,1:])]
# mesh = mymesh.Mesh(tria.points, cells)
#mesh.write("../new_mesh.stl", file_format="stl")
# mymesh.write_points_cells("../test_mesh.stl",tria.points, cells, file_format="stl")
# plt.figure
# plt.triplot(tria.points[:,0],tria.points[:,1],tria.simplices[:,1:])


# Update coordinates
utils.update_act_coords_on_ring(dm, 0)
utils.update_act_coords_on_ring(dm, 2)

# Apply flat and show actuators
dm.apply_flat()
dm.get_position(1)
dm.surface()

# # CapSens matrix
# meas_gap = utils.capsens_measure(dm, 1)