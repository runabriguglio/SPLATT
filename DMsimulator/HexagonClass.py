import numpy as np
import os
from scipy.sparse import csr_matrix
from scipy.interpolate import griddata

from read_config import readConfig
# from zernike_polynomials import computeZernike as czern
# from rotate_coordinates import cw_rotate as crot
from read_and_write_fits import write_to_fits
from read_and_write_fits import read_fits 

from rotate_coordinates import rotate_by_60deg as rot60
from rotate_coordinates import cw_rotate

# Useful variables
SIN60 = np.sin(np.pi/3.)
COS60 = np.cos(np.pi/3.)

def sum_n(n):
    return int(n*(n+1)/2)


def semi_structured_point_cloud(points_per_side):
    """ Defines a structured point cloud with random 
    displacements to act as a point mesh for the hexagon """

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


class Hexagon():
    """ Class defining the single hexagonal segment parameters """

    def __init__(self, TN):

        # Read configuration files
        config_path = './config_par_' + TN + '.yaml'
        dm_par, opt_par = readConfig(config_path)

        self.hex_len = dm_par[1]
        self.pix_scale = opt_par[0]
        self.act_pitch = dm_par[3]
        self.act_radius = dm_par[4]

        points_per_side = 9

        # Manage save path
        self.savepath = './' + TN + '/'
        try: # Create new folder
            os.mkdir(TN)
        except FileExistsError: # Folder already existing
            pass

        # Run initialization functions
        self._define_mask()
        self._initialize_act_coords()
        self._define_mesh(points_per_side)

    
    def simulate_influence_functions(self):
        """ Simulate the actuator influence functions 
        via grid interpolation """
        
        n_w = np.size(self.local_mask) # length of the masked flattened image
        n_acts = len(self.local_act_coords)
        iff_mat_shape = np.array([n_w, n_acts])
        
        file_path = self.savepath + 'fake_local_influence_function_matrix.fits'
        try:
            self.IFF = read_fits(file_path, sparse_shape = iff_mat_shape)
            return
        except FileNotFoundError:
            pass

        coords = self.local_act_coords
        x, y = coords[:, 0], coords[:, 1]
        npix_x, npix_y = np.shape(self.local_mask)  
        new_x = np.linspace(min(x), max(x), npix_y) # x coordinate is the img column!
        new_y = np.linspace(min(y), max(y), npix_x) # y coordinate is the img row!
        gx, gy = np.meshgrid(new_x, new_y)
        
        act_data = np.zeros([n_acts,n_acts])
        act_data[np.arange(n_acts),np.arange(n_acts)] = 1
        img_cube = griddata((x, y), act_data, (gx, gy), fill_value = 0.)
        
        cube_mask = np.tile(self.local_mask,n_acts)
        cube_mask = np.reshape(cube_mask, np.shape(img_cube), order = 'F')
        masked_cube = np.ma.masked_array(img_cube,cube_mask)
        
        mask = self.local_mask.copy()
        flat_mask = mask.flatten()
        ids = np.arange(len(flat_mask))
        valid_ids = ids[ids[~flat_mask]]
        row_ids = np.tile(valid_ids, n_acts)
        pix_ids = row_ids.flatten()

        act_ids = np.arange(n_acts)
        # act_ids = np.tile(act_ids, len(valid_ids))
        act_ids = np.repeat(act_ids, len(valid_ids))
        
        data = np.zeros(len(valid_ids)*n_acts)
        for i in np.arange(n_acts):
            img = img_cube[:,:,i]
            data[i*len(valid_ids):(i+1)*len(valid_ids)] = img[~self.local_mask]

        iff_mat = csr_matrix((data, (pix_ids, act_ids)), iff_mat_shape)
            
        # Save local IFF sparse mat
        self.IFF = iff_mat
        data_list = []
        data_list.append(iff_mat.data)
        data_list.append(iff_mat.indices)
        data_list.append(iff_mat.indptr)
        write_to_fits(data_list, file_path)
        
        # return masked_cube

    
    def update_act_coords(self, new_coords):
        """ Changes the local actuator coordinates 
        on the hexagonal segment to new_coords """

        self.local_act_coords = new_coords
        # The IFFs and mesh now need to be updated:
        self._simulate_influence_functions()
        self._define_mesh() 


    def _define_mask(self):
        """ Forms the hexagonal mask to be placed 
        at the given hexagonal coordinates """

        file_path = self.savepath + 'hexagon_mask.fits'
        try:
            self.local_mask = read_fits(file_path, is_bool = True)
            return
        except FileNotFoundError:
            pass

        L = self.hex_len*self.pix_scale

        # Consider points in the upper-left of the hexagon
        max_y = int(L*SIN60)
        max_x = int(L/2.+L*COS60)

        mask_ul = np.fromfunction(lambda i,j: j >= L/2. + i/SIN60*COS60, [max_y,max_x])

        mask = np.zeros([2*max_y,2*max_x], dtype=bool)

        mask[0:max_y,max_x:] = mask_ul # upper left
        mask[0:max_y,0:max_x] = np.flip(mask_ul,1) # upper right
        mask[max_y:,:] = np.flip(mask[0:max_y,:],0) # lower
                    
        self.local_mask = mask

        # Save to fits
        write_to_fits((self.local_mask).astype(np.uint8), file_path)


    def _initialize_act_coords(self):
        """ Defines the local actuator coordinates 
        on the hexagonal the segment """

        file_path = self.savepath + 'local_act_coords.fits'
        try:
            self.local_act_coords = read_fits(file_path)
            return
        except FileNotFoundError:
            pass
        
        # Normalize quantities by hexagon side length (hex_len)
        act_rad = self.act_radius/self.hex_len
        act_spacing = self.act_pitch/self.hex_len
        
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
        self.local_act_coords = act_coords*self.hex_len
        write_to_fits(act_coords, file_path)

        # self._simulate_influence_functions()
        # self._simulate_FF_matrix(CapSens_R_in = 0., 
        #                          CapSens_R_out = self.act_radius)


    def _define_mesh(self, points_per_side):
        """ Defines the local mesh point coordinates
        on the segment starting from a semi-structured
        point cloud and adding the act locations to it"""
        
        # Load or create the hexagon's mesh
        mesh_path = self.savepath + str(points_per_side) + 'pps_local_mesh_point_coords.fits'
        try:
            self.local_mesh_coords = read_fits(mesh_path)
            return
        except FileNotFoundError:
            pass

        # Load or create a 'random' mesh for the hexagon 
        file_path = self.savepath + str(points_per_side) + 'pps_point_cloud_coords.fits'
        try:
            self.point_cloud = read_fits(file_path)
        except FileNotFoundError:
            self.point_cloud = semi_structured_point_cloud(points_per_side)
            write_to_fits(self.point_cloud, file_path)
        
        # Add points on the actuator locations
        n_acts = len(self.local_act_coords)
        flat_act_coords = np.tile(self.local_act_coords,5).flatten()
        up_down_left_right = np.array([[0,0],[0,1],[0,-1],[-1,0],
                                       [1,0]]).flatten()
        UDLR = np.tile(up_down_left_right,n_acts)
        act_points = flat_act_coords + UDLR*self.act_radius/self.hex_len
        act_points = np.reshape(act_points,[n_acts*5,2])
        
        mesh_points = np.concatenate((self.point_cloud, act_points))
        
        # Save result
        self.local_mesh_coords = mesh_points
        write_to_fits(mesh_points, mesh_path)

        # self._compute_stiffness_matrix()
        # self._compute_thermal_sensitivity_matrix()

            