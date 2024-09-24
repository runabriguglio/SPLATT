import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

from rotate_coordinates import rotate_by_60deg as rot60
from rotate_coordinates import cw_rotate

def sum_n(n):
    return int(n*(n+1)/2)

# Useful variables
SIN60 = np.sin(np.pi/3.)
COS60 = np.cos(np.pi/3.)


def segment_act_coordinates(act_pitch, capsens_radius, hex_side_len):
    acts_per_side = int((hex_side_len-2*capsens_radius)/act_pitch)
    n_acts_tri = sum_n(acts_per_side)
    act_coords = np.zeros([n_acts_tri*6+1,2])
    dx = (1.-2*capsens_radius/hex_side_len)/acts_per_side
    
    for k in range(acts_per_side):
        y = np.linspace(SIN60*k*dx,0.,k+1)
        x = (k+1)*dx - COS60/SIN60 * y
        n = k+1
        p = np.zeros([6*n,2])
        p[0:n,0] = x
        p[0:n,1] = y
        p[n:2*n,:] = rot60(p[0:n,:].T)
        p[2*n:3*n,:] = rot60(p[n:2*n,:].T)
        p[3*n:,:] = cw_rotate(p[0:3*n,:].T,np.array([np.pi]))
        act_coords[1+sum_n(k)*6:1+sum_n(k+1)*6,:] = p
        
    return act_coords


def semi_structured_point_cloud(points_per_side):
    # Upper left triangle of the hexagon
    ul_triangle = np.zeros([3,2])
    ul_triangle[1,:] = np.array([2*COS60, 0.])
    ul_triangle[2,:] = np.array([0.5, SIN60])
    
    plist = np.zeros([sum_n(points_per_side+1)+3,2])

    # Structured mesh + noise
    plist[0:3,:] = ul_triangle
    dx = 1./points_per_side 
    sig = dx/7.5
    for k in np.arange(1,points_per_side+1):
        y = np.linspace(0.,SIN60*k*dx,k+1)
        x = k*dx - COS60/SIN60 * y
        if k < points_per_side:
            y[1:-1] = y[1:-1] + np.random.randn(k-1)*sig
            x[1:-1] = x[1:-1] + np.random.randn(k-1)*sig
        plist[3+sum_n(k):3+sum_n(k+1),0] = x
        plist[3+sum_n(k):3+sum_n(k+1),1] = y

    points = np.zeros([len(plist)*6,2])
    points[0:len(plist),:] = plist

    for i in range(5):
        points[(i+1)*len(plist):(i+2)*len(plist),:] = rot60(points[i*len(plist):(i+1)*len(plist),:])

    return points


def mesh_from_act_coordinates(points_per_side,act_coords,capsens_radius,hex_side_len):
    point_cloud = semi_structured_point_cloud(points_per_side)
    
    n_acts = len(act_coords)
    flat_act_coords = np.tile(act_coords,5).flatten()
    up_down_left_right = np.array([[0,0],[0,1],[0,-1],[-1,0],
                                   [1,0]]).flatten()
    UDLR = np.tile(up_down_left_right,n_acts)
    capsens_points = flat_act_coords + UDLR*capsens_radius/hex_side_len
    capsens_points = np.reshape(capsens_points,[n_acts*5,2])
    
    points = np.concatenate((point_cloud,capsens_points))
    
    return points

def distribute_local_to_global(local_coords, hex_side_len, hex_center):
    
    global_coords = local_coords * hex_side_len + hex_center
    
    return global_coords
    
   




