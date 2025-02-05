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
from scipy.spatial import Delaunay
import matrix_calculator as matcalc

coords = dm.segment[0].act_coords
pps = 40
act_radius = dm.geom.hex_side_len

points = matcalc._define_mesh(coords, pps, act_radius)
R = 5 # [m]
X_coords = points[:,0]
Y_coords = points[:,1]
Z_coords = np.sqrt(R**2 - (X_coords**2 + Y_coords**2))
points = np.vstack((X_coords, Y_coords, Z_coords))
points = np.transpose(points, (1,0))
tria = Delaunay(points)

import meshio as mymesh

cells = [("triangle", tria.simplices[:,1:])]

# mesh = meshio.Mesh(points, cells)
# mesh.write("../mesh.vtk", file_format="vtk")

mymesh.write_points_cells("../test_mesh.stl",tria.points, cells, file_format="stl")


# Update coordinates
utils.update_act_coords_on_ring(dm, 0)
utils.update_act_coords_on_ring(dm, 2)

# Apply flat and show actuators
dm.apply_flat()
dm.get_position(1)
dm.surface()

# # CapSens matrix
# meas_gap = utils.capsens_measure(dm, 1)