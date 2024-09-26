import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

from read_config import readConfig
from read_and_write_fits import write_to_fits
from read_and_write_fits import read_fits_file as read_fits

from rotate_coordinates import rotate_by_60deg as rot60
from rotate_coordinates import cw_rotate

def sum_n(n):
    return int(n*(n+1)/2)

# Useful variables
SIN60 = np.sin(np.pi/3.)
COS60 = np.cos(np.pi/3.)


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

    


class Mesh():

    def __init__(self, TN):

        # Read configuration files
        config_path = './config_par_' + TN + '.yaml'
        dm_par, opt_par = readConfig(config_path)

        self.hex_side_len = dm_par[1]
        self.act_pitch = dm_par[3]
        self.act_radius = dm_par[4]

        self.pix_scale = opt_par[0]
        
        self.savepath = './' + TN + '/'
        
        # try: # Create new folder
        #     os.mkdir(TN)
        # except FileExistsError: # Folder already existing
        #     pass
            
        self.segment_act_coordinates()
        self.mesh_from_act_coordinates(9)
        # tria = Delaunay(mesh_points)
        
        
        
    def segment_act_coordinates(self):
        
        file_path = self.savepath + 'local_act_coords.fits'
        try:
            self.local_act_coords = read_fits(file_path)
            return
        
        except FileNotFoundError:
            pass
        
        # Normalize quantities by hexagon side
        act_rad = self.act_radius/self.hex_side_len
        act_spacing = self.act_pitch/self.hex_side_len
        
        acts_per_side = int((1.-2*act_rad)/act_spacing)
        n_acts_tri = sum_n(acts_per_side)
        
        act_coords = np.zeros([n_acts_tri*6+1,2])
        dx = (1.-2*act_rad)/acts_per_side
        
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
            
        # Save result
        self.local_act_coords = act_coords
        write_to_fits(act_coords, file_path)
        
        
        
    def mesh_from_act_coordinates(self, points_per_side):
        
        file_path = self.savepath + 'local_mesh_points_coords.fits'
        try:
            self.local_mesh_coords = read_fits(file_path)
            return
        
        except FileNotFoundError:
            pass
        
        # Generate a semi-structured mesh
        point_cloud = semi_structured_point_cloud(points_per_side)
        
        # Add points on the actuator locations
        n_acts = len(self.local_act_coords)
        flat_act_coords = np.tile(self.local_act_coords,5).flatten()
        up_down_left_right = np.array([[0,0],[0,1],[0,-1],[-1,0],
                                       [1,0]]).flatten()
        UDLR = np.tile(up_down_left_right,n_acts)
        act_points = flat_act_coords + UDLR*self.act_radius/self.hex_side_len
        act_points = np.reshape(act_points,[n_acts*5,2])
        
        mesh_points = np.concatenate((point_cloud,act_points))
        
        # Save result
        self.local_mesh_coords = mesh_points
        write_to_fits(mesh_points, file_path)
        
        return mesh_points
    
    
    
    def act_influence_functions(self, mask):
        
        file_path = self.savepath + 'fake_local_influence_functions.fits'
        try:
            self.local_IFFs = read_fits(file_path, list_len = len(self.local_act_coords), is_ma = True)
            return
        
        except FileNotFoundError:
            pass
        
        dist = lambda x,y: np.sqrt(x**2+y**2)

        X,Y = np.shape(mask)
        r = np.max([X,Y])
        
        masked_img_cube = []
        hex_pix_len = self.hex_side_len * self.pix_scale
        
        for k, coords in enumerate(self.local_act_coords):
            x_act = int(coords[0]*hex_pix_len) + Y/2
            y_act = int(coords[1]*hex_pix_len) + X/2
            
            img = np.fromfunction(lambda i,j: r - dist(j-x_act, i-y_act), [X,Y])
            
            # Normalization to uint8
            masked_data = img[~mask]
            norm = np.max(masked_data)-np.min(masked_data)
            uint8_img = ((2**8-1) * (img - np.min(masked_data))/norm).astype(np.uint8)
            
            masked_img = np.ma.masked_array(uint8_img, mask)
            masked_img_cube.append(masked_img)
            
            # masked_data = masked_img.data[~masked_img.mask]
            # print([np.max(masked_data),np.min(masked_data)])
            
            # plt.figure()
            # plt.imshow(masked_img,origin='lower')
            # plt.scatter(x_act,y_act,c='black')
            
        # Save result
        self.local_IFFs = masked_img_cube
        write_to_fits(masked_img_cube, file_path)
            
        return masked_img_cube





