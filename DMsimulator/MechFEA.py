
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import numpy as np
from finite_elements import semi_structured_point_cloud


mesh_points = semi_structured_point_cloud(80)
tria = Delaunay(mesh_points)

# Material properties
E = 90.3e9 #[Pa]
nu = 0.24
a = E/(1-nu**2)
b = nu*a
c = E/2/(1+nu)

E_mat = np.array([[a,b,0],[b,a,0],[0,0,c]])

# Shell's thickness
t = 1e-3 #[m]

stiff_mat = t * E_mat
bend_mat = t**3/12 * E_mat