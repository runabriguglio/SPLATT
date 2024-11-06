import numpy as np
import os
from scipy.sparse import csr_matrix
# from scipy.interpolate import griddata

from read_config import readConfig
from zernike_polynomials import computeZernike as czern
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

        self.points_per_side = 9

        # Manage save path
        self.savepath = './' + TN + '/'
        try: # Create new folder
            os.mkdir(TN)
        except FileExistsError: # Folder already existing
            pass

        # Run initialization functions
        self._define_mask()
        self._initialize_act_coords()
        self._define_mesh()


    def compute_interaction_matrix(self, n_modes):
        """ Computes the interaction matrix: 
            [n_pixels,n_hexes*n_modes] """

        int_mat_shape = [np.size(self.local_mask),n_modes]
            
        file_path = self.savepath + str(n_modes) + 'modes_local_interaction_matrix.fits'
        try:
            self.int_mat = read_fits(file_path, sparse_shape = int_mat_shape)
            return
        except FileNotFoundError:
            pass
        
        hex_data_len = np.sum(1-self.local_mask)
        modes = np.zeros([hex_data_len*n_modes])
        for j in range(n_modes):
            modes[hex_data_len*j:hex_data_len*(j+1)] = czern(j+1, self.local_mask)
            
        print('Computing interaction matrix...')      
        row_indices = np.tile(self.hex_valid_ids, n_modes)
        row = row_indices.flatten()
        
        mode_indices = np.arange(n_modes)
        col = np.repeat(mode_indices, hex_data_len)

        self.int_mat = csr_matrix((modes, (row,col)),  
                                int_mat_shape)
        
        # Save to fits
        print('Saving interaction matrix...') 
        data_list = []
        data_list.append(self.int_mat.data)
        data_list.append(self.int_mat.indices)
        data_list.append(self.int_mat.indptr)
        write_to_fits(data_list, file_path)

    
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

        # Save valid indices
        flat_mask = mask.flatten()
        ids = np.arange(np.size(mask))
        self.hex_valid_ids = ids[ids[~flat_mask]]

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

        # Compute influence functions
        self._simulate_influence_functions()

        # self._simulate_FF_matrix(CapSens_R_in = 0., 
        #                          CapSens_R_out = self.act_radius)


    def _define_mesh(self):
        """ Defines the local mesh point coordinates
        on the segment starting from a semi-structured
        point cloud and adding the act locations to it"""

        # Load or create a 'random' mesh for the hexagon 
        file_path = self.savepath + 'point_cloud_coords.fits'
        try:
            self.point_cloud = read_fits(file_path)
        except FileNotFoundError:
            point_cloud = semi_structured_point_cloud(self.points_per_side)
            write_to_fits(point_cloud, file_path)

        # Load or create the hexagon's mesh
        mesh_path = self.savepath + 'local_mesh_points_coords.fits'
        try:
            self.local_mesh_coords = read_fits(mesh_path)
            return
        except FileNotFoundError:
            pass
        
        # Add points on the actuator locations
        n_acts = len(self.local_act_coords)
        flat_act_coords = np.tile(self.local_act_coords,5).flatten()
        up_down_left_right = np.array([[0,0],[0,1],[0,-1],[-1,0],
                                       [1,0]]).flatten()
        UDLR = np.tile(up_down_left_right,n_acts)
        act_points = flat_act_coords + UDLR*self.act_radius/self.hex_len
        act_points = np.reshape(act_points,[n_acts*5,2])
        
        mesh_points = np.concatenate((self.point_cloud,act_points))
        
        # Save result
        self.local_mesh_coords = mesh_points
        write_to_fits(mesh_points, mesh_path)

        # self._compute_stiffness_matrix
        # self._compute_thermal_sensitivity_matrix


    def _simulate_influence_functions(self):
        """ Simulate the actuator influence functions 
        via grid interpolation """
        
        n_w = len(self.hex_valid_ids) # length of the masked flattened image
        n_acts = len(self.local_act_coords)
        iff_mat_shape = np.array([n_w, n_acts])
        
        file_path = self.savepath + 'fake_local_influence_function_matrix.fits'
        try:
            self.IFF = read_fits(file_path, sparse_shape = iff_mat_shape)
            return
        except FileNotFoundError:
            pass

        # print('Computing IFF matrix...')      
        # #x, y = in_mesh[:, 0], in_mesh[:, 1]
        # #new_x = np.linspace(min(x), max(x), npix)
        # #new_y = np.linspace(min(y), max(y), npix)
        # #gx, gy = np.meshgrid(new_x, new_y)
        # #if_grid = griddata((x, y), if_matteo, (gx, gy), method='linear'

        # act_data = np.zeros(n_acts*n_acts)
        # act_data[n_acts*np.arange(n_acts)] = 1


        # row_ids = np.tile(self.hex_valid_ids, n_acts)
        # pix_ids = row_ids.flatten()

        # act_ids = np.arange(n_acts)

        # iff_mat = csr_matrix((data, (pix_ids, act_ids)), iff_mat_shape)
            
        # # Save local IFF sparse mat
        # print('Saving influence function matrix...') 
        # self.local_IFF = iff_mat
        # data_list = []
        # data_list.append(iff_mat.data)
        # data_list.append(iff_mat.indices)
        # data_list.append(iff_mat.indptr)
        # write_to_fits(data_list, file_path)
            