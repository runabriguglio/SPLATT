import numpy as np
import os
from scipy.interpolate import griddata
from tps import ThinPlateSpline # for the simulated IFF

from read_configuration import read_config
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


class Calculator():
    """ Class defining the single hexagonal segment parameters """

    def __init__(self, TN):

        # Read configuration files
        config_path = './config_par_' + TN + '.yaml'
        dm_par, opt_par = read_config(config_path)

        self.hex_len = dm_par[1]
        self.pix_scale = opt_par[0]
        self.act_pitch = dm_par[3]
        self.act_radius = dm_par[4]

        # Manage save path
        self.savepath = './' + TN + '/'
        try: # Create new folder
            os.mkdir(TN)
        except FileExistsError: # Folder already existing
            pass

        # Run initialization functions
        self._define_mask()
        self._initialize_act_coords()
        
    @staticmethod
    def compute_reconstructor(IFF):
        return np.linalg.pinv(IFF)
    
    def compute_influence_functions(self, simulate_IFFs = False, act_coords = None):
        """ Wrapper to simulate/compute the influence functions """
        
        if simulate_IFFs:
            self.simulate_influence_functions(act_coords)
        else:
            self.simulate_influence_functions(act_coords)
        
    
    def calculate_influence_functions(self, act_coords = None):
        """ Project the actuator influence functions 
        on the mask via grid interpolation """
        
        file_path = self.savepath + 'local_influence_function_matrix.fits'
        try:
            self.IFF = read_fits(file_path)
        except FileNotFoundError:
            pass

        if act_coords is None:
            act_coords = self._initialize_act_coords()
        n_acts = len(act_coords)
        
        # IFF data on mesh
        points_per_side = 9
        self._define_mesh(act_coords, points_per_side)
        # iff_data = self._compute_iff_data()
        
        # Grid interpolation
        x, y = act_coords[:, 0], act_coords[:, 1]
        npix_x, npix_y = np.shape(self.local_mask)  
        new_x = np.linspace(min(x), max(x), npix_y) # x coordinate is the img column!
        new_y = np.linspace(min(y), max(y), npix_x) # y coordinate is the img row!
        gx, gy = np.meshgrid(new_x, new_y)
        img_cube = griddata((x, y), self.iff_data, (gx, gy), fill_value = 0., method = 'linear')
        
        # Masked array
        cube_mask = np.tile(self.local_mask,n_acts)
        cube_mask = np.reshape(cube_mask, np.shape(img_cube), order = 'F')
        cube = np.ma.masked_array(img_cube,cube_mask)
        
        # Save image cube to fits
        write_to_fits(cube, self.savepath + 'influence_functions_image_cube.fits')
        
        # Save valid data to IFF (full) matrix
        flat_cube = cube.data[~cube.mask]
        valid_len = np.sum(1-self.local_mask)
        local_IFF = np.reshape(flat_cube, [valid_len, n_acts])
        
        self.IFF = np.array(local_IFF)
        
        # Save to fits
        write_to_fits(self.IFF, file_path)
        
        
    def simulate_influence_functions(self, act_coords = None):
        """ Simulate the influence functions by 
        imposing 'perfect' zonal commands """
        
        file_path = self.savepath + 'local_simulated_influence_function_matrix.fits'
        try:
            self.sim_IFF = read_fits(file_path)
            return
        except FileNotFoundError:
            pass
        
        if act_coords is None:
            act_coords = self._initialize_act_coords()
        
        n_acts = len(act_coords)
        valid_len = np.sum(1-self.local_mask)
        
        max_x, max_y = np.shape(self.local_mask)
        
        pix_coords = np.zeros([max_x*max_y,2])
        pix_coords[:,0] = np.repeat(np.arange(max_x),max_y)
        pix_coords[:,1] = np.tile(np.arange(max_y),max_x)
        
        act_pix_coords = np.zeros([n_acts,2])
        act_pix_coords[:,0] = (act_coords[:,1] * self.pix_scale + max_x/2).astype(int)
        act_pix_coords[:,1] = (act_coords[:,0] * self.pix_scale + max_y/2).astype(int)
        
        img_cube = np.zeros([max_x,max_y,n_acts])

        for k in range(n_acts):
            act_data = np.zeros(n_acts)
            act_data[k] = 1
            tps = ThinPlateSpline(alpha=0.0)
            tps.fit(act_pix_coords, act_data)
            flat_img = tps.transform(pix_coords)
            img_cube[:,:,k] = np.reshape(flat_img, [max_x,max_y])

        # Masked array
        cube_mask = np.tile(self.local_mask,n_acts)
        cube_mask = np.reshape(cube_mask, np.shape(img_cube), order = 'F')
        cube = np.ma.masked_array(img_cube,cube_mask)
        
        # Save image cube to fits
        write_to_fits(cube, self.savepath + 'simulated_influence_functions_image_cube.fits')

        # Save valid data to IFF (full) matrix
        flat_cube = cube.data[~cube.mask]
        local_IFF = np.reshape(flat_cube, [valid_len, n_acts])
        
        self.sim_IFF = np.array(local_IFF)
        
        # Save to fits
        write_to_fits(self.sim_IFF, file_path)
        
        
    def _define_mask(self):
        """ Forms the hexagonal mask to be placed 
        at the given hexagonal coordinates """

        file_path = self.savepath + 'local_mask.fits'
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
        write_to_fits((mask).astype(np.uint8), file_path)


    def _initialize_act_coords(self):
        """ Defines the local actuator coordinates 
        on the hexagonal the segment """

        file_path = self.savepath + 'local_act_coords.fits'
        try:
            self.local_act_coords = read_fits(file_path)
            return self.local_act_coords
        except FileNotFoundError:
            pass
        
        # Normalize quantities by hexagon side length (hex_len)
        L = self.hex_len
        rad = self.act_radius/L
        pitch = self.act_pitch/L
        
        acts_per_side = (1+pitch)/(2*rad + pitch) #+ 1
        dx = 2*rad+pitch
        
        acts_per_side = int(acts_per_side)
        n_acts_tri = sum_n(acts_per_side)
        
        act_coords = np.zeros([n_acts_tri*6+1,2])
        
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
            
        # Rescaling
        self.local_act_coords = act_coords*self.hex_len #*(acts_per_side)/(acts_per_side-1)
            
        # Save result
        write_to_fits(self.local_act_coords, file_path)
        
        return self.local_act_coords


    def _define_mesh(self, act_coords, points_per_side):
        """ Defines the local mesh point coordinates
        on the segment starting from a semi-structured
        point cloud and adding the act locations to it"""
        
        # Load or create the hexagon's mesh
        mesh_path = self.savepath + 'local_mesh_point_coords.fits'

        # Load or create a 'random' mesh for the hexagon 
        file_path = self.savepath + str(points_per_side) + 'pps_point_cloud_coords.fits'
        try:
            self.point_cloud = read_fits(file_path)
        except FileNotFoundError:
            self.point_cloud = semi_structured_point_cloud(points_per_side)
            write_to_fits(self.point_cloud, file_path)
        
        # Add points on the actuator locations
        n_acts = len(act_coords)
        flat_act_coords = np.tile(act_coords,5).flatten()
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
        # self._compute_iff_data
        # self._compute_thermal_sensitivity_matrix()
        
        

            