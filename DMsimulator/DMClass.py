import numpy as np
import os
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt

from read_config import readConfig
from zernike_polynomials import computeZernike as czern
from rotate_coordinates import cw_rotate as crot
from read_and_write_fits import write_to_fits
from read_and_write_fits import read_fits 

# Useful variables
SIN60 = np.sin(np.pi/3.)
COS60 = np.cos(np.pi/3.)

def n_hexagons(n_rings):
    return int(1 + (6 + n_rings*6)*n_rings/2)


class Segment():
    """ Class defining the single segment
    with center coordinates and actuator displacements """
    
    def __init__(self, segment_id, center_coordinates, segment_mask):
        
        self.id = segment_id
        self.mask = segment_mask
        self.center = center_coordinates
        self.wavefront = np.zeros(np.sum(1-self.mask))
        
        self.act_pos = []
        
        
    def wavefront(self):
        """ Plots an image of the segment's current shape 
        on the local segment mask """
        
        # Project wavefront on mask
        mask = self.mask.copy()
        flat_mask = mask.flatten()
        image = np.zeros(np.size(flat_mask))
        image[~flat_mask] = self.wavefront
        image = np.reshape(image, np.shape(mask))
        image = np.ma.masked_array(image, mask)
        
        # Plot image
        plt.figure()
        plt.imshow(image, origin = 'lower', cmap = 'inferno')
        plt.colorbar()
        plt.title('Segment ' + str(self.id) + ' shape')
        
        
    def get_position(self):
        """ Wrapper to read the current 
        position of the segments's actuators """
        
        pos = self.act_pos
        return pos
    
    
    def set_position(self, pos_cmd, absolute_pos = True):
        """ Command the position of the 
        segments's actuators in absolute (default) or relative terms"""
        
        old_pos = self.act_pos
        
        if absolute_pos:
            new_pos = pos_cmd
        else:
            new_pos = pos_cmd + old_pos

        self.act_pos = new_pos
        # self.wavefront += IFF @ (new_pos-old_pos)
        
        return new_pos
    
    
    # def flat_cmd(self):
        
    #     # Compute cmd and wavefront change
    #     act_cmd = R @ self.wavefront
    #     delta_wf = IM @ act_cmd
        
    #     # Update act position and wavefront
    #     self.act_pos += act_cmd
    #     self.wavefront += delta_wf
        
    #     flat_rms = np.std(self.wavefront)
        
    #     return flat_rms#, act_cmd
        
        
        

class DM():
    """ Class defining the segmented deformable mirror """

    def __init__(self, TN):
        
        # Read configuration files
        config_path = './config_par_' + TN + '.yaml'
        dm_par, opt_par = readConfig(config_path)
        
        self.gap = dm_par[0]
        self.hex_side_len = dm_par[1]
        self.n_rings = int(dm_par[2])

        self.pix_scale = opt_par[0]
        
        self.savepath = './' + TN + '/'
        self.n_hex = n_hexagons(self.n_rings)
        
        try: # Create new folder
            os.mkdir(TN)
        except FileExistsError: # Folder already existing
            pass
        
        self._compute_segment_centers()
        self._define_global_valid_ids()
        self._assemble_global_mask()
        self._define_segment_array()
        
        
    def plot_wavefront(self, wavefront, title_str = None):
        """ Plots an image on the global mask 
        given a wavefront """
        
        full_img = np.reshape(wavefront, np.shape(self.global_mask))
        img = np.ma.masked_array(full_img, self.global_mask)

        plt.figure()
        plt.imshow(img, origin = 'lower', cmap='inferno')
        if title_str is not None:
            plt.title(title_str)
            
        # # Add a colorbar
        # plt.colorbar()
        # fig.set_clim([min(wavefront), max(wavefront)])
        
        
    def segment_scramble(self):
        """ Defines an initial segment scramble
        returning a masked array image """
        
        file_path = self.savepath + 'segment_scramble.fits'
        try:
            masked_img = read_fits(file_path, is_ma = True)
            return masked_img
        except FileNotFoundError:
            pass
        
        n_hex = n_hexagons(self.n_rings)
        
        # Retrieve number of modes from the interaction matrix
        n_modes = int(np.shape(self.int_mat)[1]/n_hex)
        
        # Generate random mode coefficients
        mode_vec = np.random.randn(n_hex*n_modes)
        
        # Probability inversely proportional to spatial frequency
        m = int(np.ceil((np.sqrt(8*n_modes)-1.)/2.))
        freq_vec = np.repeat(np.arange(m)+1,np.arange(m)+1)
        prob_vec = 1./freq_vec[:n_modes]
        prob_vec_rep = np.tile(prob_vec,n_hex)
        
        # Modulate on the probability
        mode_vec = mode_vec * prob_vec_rep
        
        # Matrix product
        flat_img = self.int_mat*mode_vec
        
        # # Global modes
        # n_glob_modes = np.shape(self.glob_int_mat)[1]
        # glob_mode_vec = np.random.randn(n_glob_modes)
        # flat_img += self.glob_int_mat*glob_mode_vec
        
        # Reshape and mask image
        img = np.reshape(flat_img, np.shape(self.global_mask))
        masked_img = np.ma.masked_array(img, self.global_mask)
        
        # Save to fits
        write_to_fits(masked_img, file_path)
        return masked_img
        
        
    def compute_global_interaction_matrix(self,n_modes):
        """ Computes the global interaction matrix: 
            [n_pixels,n_modes] """
            
        int_mat_shape = [np.size(self.global_mask),n_modes]
            
        file_path = self.savepath + str(n_modes) + 'modes_global_interaction_matrix.fits'
        try:
            self.glob_int_mat = read_fits(file_path, sparse_shape = int_mat_shape)
            return
        except FileNotFoundError:
            pass
        
        data_len = np.sum(1-self.global_mask)
        modes_data = np.zeros([data_len*n_modes])
        for j in range(n_modes):
            modes_data[data_len*j:data_len*(j+1)] = czern(j+1, self.global_mask)
            
        mask = self.global_mask.copy()
        flat_mask = mask.flatten()
        ids = np.arange(len(flat_mask))
        valid_ids = ids[ids[~flat_mask]]
        row_indices = np.tile(valid_ids,n_modes)
        row = row_indices.flatten()
        
        mode_indices = np.arange(n_modes)
        col = np.repeat(mode_indices,data_len)
    
        self.glob_int_mat = csr_matrix((modes_data, (row,col)),  
                                  int_mat_shape)
        
        # Save to fits
        data_list = []
        data_list.append(self.glob_int_mat.data)
        data_list.append(self.glob_int_mat.indices)
        data_list.append(self.glob_int_mat.indptr)
        write_to_fits(data_list, file_path)
        
        
        
    def compute_interaction_matrix(self,n_modes):
        """ Computes the interaction matrix: 
            [n_pixels,n_hexes*n_modes] """
            
        n_hex = n_hexagons(self.n_rings)
        int_mat_shape = [np.size(self.global_mask),n_modes*n_hex]
            
        file_path = self.savepath + str(n_modes) + 'modes_interaction_matrix.fits'
        try:
            self.int_mat = read_fits(file_path, sparse_shape = int_mat_shape)
            return
        except FileNotFoundError:
            pass
        
        hex_data_len = np.sum(1-self.local_mask)
        row_modes = np.zeros([hex_data_len*n_modes])
        for j in range(n_modes):
            row_modes[hex_data_len*j:hex_data_len*(j+1)] = czern(j+1, self.local_mask)
            
        row_indices = np.tile(self.hex_valid_ids,n_modes)
        row = row_indices.flatten()
        
        mode_indices = np.arange(int(n_modes*n_hex))
        col = np.repeat(mode_indices,hex_data_len)
        
        data = np.tile(row_modes,n_hex)
    
        self.int_mat = csr_matrix((data, (row,col)),  
                                  int_mat_shape)
        
        # Save to fits
        data_list = []
        data_list.append(self.int_mat.data)
        data_list.append(self.int_mat.indices)
        data_list.append(self.int_mat.indptr)
        write_to_fits(data_list, file_path)
        
        
    def _assemble_global_mask(self):
        """ Assembles the global segmented mask from
        the local mask data """
        
        file_path = self.savepath + 'global_mask.fits'
        try:
            self.global_mask = read_fits(file_path, is_bool = True)
            return
        except FileNotFoundError:
            pass
        
        # Height of hex + gap
        L = self.gap + 2.*self.hex_side_len*SIN60
        
        # Full mask dimensions
        Ny = (L*self.n_rings + self.hex_side_len*SIN60)*2*self.pix_scale
        Nx = (L*self.n_rings*SIN60 + self.hex_side_len*(0.5+COS60))*2*self.pix_scale
        mask_shape = np.array([int(Ny),int(Nx)])
        mask = np.zeros(mask_shape, dtype = bool)
        
        # Hexagon centers pixel coordinates
        pix_coords = self.hex_centers*self.pix_scale + np.array([Nx,Ny])/2.
        
        valid_len = np.sum(1-self.local_mask)
        rep_pix_coords = np.repeat(pix_coords, valid_len, axis = 0)
        
        # Data
        data = np.ones(len(rep_pix_coords), dtype=bool)
        
        # Sparse matrix assembly
        sparse_mat = csr_matrix((data, (self.global_row_idx, self.global_col_idx)),  
                                  mask_shape, dtype=bool)
        mask += sparse_mat
            
        self.global_mask = np.array(~mask)
        
        # Save to fits
        write_to_fits((self.global_mask).astype(np.uint8), file_path)
        
        
        
    def _define_segment_array(self):
        """ Builds an array of Segment class objects,
        containing their center coordinates and the
        actuator positions """
        
        self.segments = []
        
        for k,coords in enumerate(self.hex_centers):
            self.segments.append(Segment(k, coords, self.local_mask))
            # local_INTMAT = self.int_mat[:,N_modes*k:N_modes*(k+1)]
            
            
    def _define_global_valid_ids(self):
        """ Finds the full aperture image (containing all segments)
        row and column indices for the segments images """
                  
        # Read local mask file
        local_mask_path = self.savepath + 'hexagon_mask.fits'
        self.local_mask = read_fits(local_mask_path, is_bool = True)
        
        file_path = self.savepath + 'segments_indices.fits'
        try:
            out = read_fits(file_path, list_len = 3)
            self.hex_valid_ids = out[0]
            self.global_row_idx = out[1]
            self.global_col_idx = out[2]
            return
        except FileNotFoundError:
            pass
        
        # Height of hex + gap
        L = self.gap + 2.*self.hex_side_len*SIN60
        
        # Full mask dimensions
        Ny = (L*self.n_rings + self.hex_side_len*SIN60)*2*self.pix_scale
        Nx = (L*self.n_rings*SIN60 + self.hex_side_len*(0.5+COS60))*2*self.pix_scale
        Ntot = np.array([Nx,Ny])

        # Hexagon centers pixel coordinates
        self.pix_coords = self.hex_centers*self.pix_scale + Ntot/2.
        
        My,Mx = np.shape(self.local_mask)
        x = np.arange(Mx,dtype=int)
        y = np.arange(My,dtype=int)
        X,Y = np.meshgrid(x,y)
        local_X = X[~self.local_mask]
        local_Y = Y[~self.local_mask]
        local_row_idx = local_Y - int(My/2)
        local_col_idx = local_X - int(Mx/2)

        n_hex = n_hexagons(self.n_rings)
        rep_local_row = np.tile(local_row_idx,n_hex)
        rep_local_col = np.tile(local_col_idx,n_hex)
        
        valid_len = np.sum(1-self.local_mask)
        rep_pix_coords = np.repeat(self.pix_coords, valid_len, axis = 0)
        
        self.global_row_idx = (rep_local_row + rep_pix_coords[:,1]).astype(int)
        self.global_col_idx = (rep_local_col + rep_pix_coords[:,0]).astype(int)
        
        # Save valid hexagon indices
        hex_idx = self.global_col_idx + self.global_row_idx*int(Nx)
        self.hex_valid_ids = np.reshape(hex_idx,[n_hex,valid_len])
        
        # Save to fits
        write_to_fits([self.hex_valid_ids,self.global_row_idx,self.global_col_idx], file_path)
        
        
    def _compute_segment_centers(self):
        """ Defines and saves the coordinates of the 
        centers of all hexagonal segments """
        
        file_path = self.savepath + 'hex_centers_coords.fits'
        try:
            self.hex_centers = read_fits(file_path)
            return
        except FileNotFoundError:
            pass
        
        # Number of hexes
        self.hex_centers = np.zeros([n_hexagons(self.n_rings),2])
        
        # Height of hex + gap
        L = self.gap + 2.*self.hex_side_len*SIN60
        angles = np.pi/3. * (np.arange(5)+1)
        
        for ring_ctr in range(self.n_rings):
            
            hex_ctr = n_hexagons(ring_ctr)
            ring_ctr += 1
            R = L*ring_ctr
            
            aux = np.array([R*SIN60,R*COS60])
            self.hex_centers[hex_ctr,:] = aux
            
            self.hex_centers[hex_ctr+ring_ctr:hex_ctr+6*ring_ctr:ring_ctr,:] = crot(aux, angles)
            
            if ring_ctr > 1:
                for j in range(ring_ctr-1):
                    shift = self.gap + 2.*self.hex_side_len*SIN60
                    aux[0] = self.hex_centers[hex_ctr,0] 
                    aux[1] = self.hex_centers[hex_ctr,1] - (j+1)*shift
                    self.hex_centers[hex_ctr+j+1,:] = aux
                    self.hex_centers[hex_ctr+j+1+ring_ctr:hex_ctr+j+1+6*ring_ctr:ring_ctr,:] = crot(aux, angles)

        # Save to fits
        write_to_fits(self.hex_centers, file_path)
        
        
        # def draw_hex_outline(self):
        #     """ Plots the hexagons' outline and the inscribed circle """
        
        #     hex_sides = np.zeros([8,2])
        #     hex_sides[0,:] = np.array([-2*COS60, 0.])
        #     hex_sides[1,:] = np.array([-0.5, SIN60])
        #     hex_sides[2,:] = np.array([-hex_sides[1,0],hex_sides[1,1]])
        #     hex_sides[3,:] = -hex_sides[0,:]
        #     hex_sides[4,:] = -hex_sides[1,:]
        #     hex_sides[5,:] = -hex_sides[2,:]
        #     hex_sides[6,:] = hex_sides[0,:]
        #     hex_sides[-1,:] = np.array([None,None])
            
        #     plt.figure()
        #     plt.grid('on')
            
        #     rep_c_coords = np.tile(self.hex_centers,len(hex_sides))
        #     rep_c_coords = rep_c_coords.flatten()
        #     hex_sides = hex_sides.flatten()
        #     rep_hex_sides = np.tile(hex_sides,len(self.hex_centers))
        #     coords = rep_c_coords + rep_hex_sides 
        #     coords = np.reshape(coords,[int(len(coords)/2),2])
        #     plt.plot(coords[:,0],coords[:,1],color='goldenrod')
                     
            
        #     # Plot inscribed cirle
        #     L = self.gap + 2.*self.hex_side_len*SIN60
        #     R = np.sqrt((L*self.n_rings)**2 + (self.hex_side_len*(0.5+COS60))**2) - self.hex_side_len*COS60
        #     x_vec = np.linspace(-R,R,100)
        #     y_vec = np.sqrt(R**2-x_vec**2)
            
        #     plt.plot(x_vec,y_vec,'--',color='green')
        #     plt.plot(x_vec,-y_vec,'--',color='green')
                
        #     plt.axis('equal')
            
        #     return coords
        