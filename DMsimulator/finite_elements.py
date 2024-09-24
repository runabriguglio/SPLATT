import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

from rotate_coordinates import rotate_by_60deg as rot60
from rotate_coordinates import cw_rotate

##########
SIN60 = np.sin(np.pi/3.)
COS60 = np.cos(np.pi/3.)

def sum_n(n):
    return int(n*(n+1)/2)



ul_triangle = np.zeros([3,2])
ul_triangle[1,:] = np.array([2*COS60, 0.])
ul_triangle[2,:] = np.array([0.5, SIN60])

# # Structured mesh
# points_per_side = 10
# plist = np.zeros([sum_n(points_per_side+1)+3,2])

# plist[0:3,:] = ul_triangle
# dx = 1./points_per_side
# for k in np.arange(1,points_per_side+1):
#     y = np.linspace(0.,SIN60*k*dx,k+1)
#     x = k*dx - COS60/SIN60 * y
#     plist[3+sum_n(k):3+sum_n(k+1),0] = x
#     plist[3+sum_n(k):3+sum_n(k+1),1] = y
    
# # Irregular mesh
# Nel = 200
# mu_x = COS60
# # mu_y = SIN60
# sigma_x = mu_x/2.4
# # sigma_y = mu_y/3
# yvec = np.random.random(Nel)
# xvec = mu_x + sigma_x * np.random.randn(Nel)

# xvec[xvec < 0.] = 1. - xvec[xvec < 0.]
# xvec[xvec > 1.] = xvec[xvec > 1.] - 1.

# xvec = np.sort(xvec)

# y_start_max = xvec[xvec <= 0.5]*SIN60/COS60
# y_end_max = (1.-xvec[xvec > 0.5])*SIN60/COS60

# y_max = np.concatenate((y_start_max,y_end_max))
# yvec[yvec > y_max] = y_max[yvec > y_max]


# plist = np.zeros([len(yvec)+3,2])
# plist[0:3,:] = ul_triangle
# plist[3:,0] = xvec
# plist[3:,1] = yvec


# # Semi-structured mesh
# points_per_side = 12
# plist = np.zeros([sum_n(points_per_side+2),2])
# ang = np.pi/3.

# plist[0:3,:] = ul_triangle
# dx = 1./points_per_side
# ctr = 3
# for k in np.arange(points_per_side):
#     p = np.zeros([k+2,2])
#     p[0,:] = [k*dx,0]
#     p[-1,:] = [k*dx*COS60,k*dx*SIN60]
#     angle = np.random.uniform(ang/(k+5),ang*(k+5)/(k+6),k) #np.random.randn(k)*sig + mu
#     angle = np.sort(angle)
#     p[1:-1,:] = cw_rotate(p[0,:],-angle)*SIN60
#     # print(3+sum_n(k)+k+2-len(plist))
#     plist[ctr:ctr+len(p),:] = p
#     ctr += len(p)
    
# y = np.linspace(0.,SIN60,k+1)
# x = 1. - COS60/SIN60 * y
# plist[ctr:,0] = x
# plist[ctr:,1] = y

# Structured mesh + noise
points_per_side = 16
plist = np.zeros([sum_n(points_per_side+1)*2+3,2])

plist[0:3,:] = ul_triangle
dx = 1./points_per_side 
sig = dx/8.
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



tria = Delaunay(points) #,incremental=True
plt.figure()
plt.grid('on')
plt.axis('equal')
plt.triplot(points[:,0], points[:,1], tria.simplices)
plt.plot(plist[:,0], plist[:,1], 'o')
plt.show()

# tria.add_points(np.array([[0.222,0.234],[0.5,0.717],[0.333,0.125],[0.6,0.1222]]),restart=True)
# plt.figure()
# plt.grid('on')
# plt.axis('equal')
# plt.triplot(points[:,0], points[:,1], tria.simplices)


