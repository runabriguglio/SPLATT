import numpy as np
import os
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt

from read_configuration import read_config
from zernike_polynomials import compute_zernike as czern
from rotate_coordinates import cw_rotate as crot
from read_and_write_fits import write_to_fits
from read_and_write_fits import read_fits 
from read_and_write_fits import write_csr_to_fits

from mirror_segments import Segment
from local_influence_functions_calculator import Calculator

# Useful variables
SIN60 = np.sin(np.pi/3.)
COS60 = np.cos(np.pi/3.)

def n_hexagons(n_rings):
    return int(1 + (6 + n_rings*6)*n_rings/2)

class SegmentedMirror():
    """ Class defining the segmented deformable mirror """

    def __init__(self, TN):
        
        # Read configuration files
        config_path = './config_par_' + TN + '.yaml'
        dm_par, opt_par = read_config(config_path)
        
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
        
        self.LMC = Calculator(TN)
        self.local_mask = self.LMC.local_mask
        
        self._compute_segment_centers()
        self._initialize_global_actuator_coordinates()
        # self._define_global_valid_ids()
        self._assemble_global_mask()
        self._define_segment_array()
        
        
    def plot_wavefront(self, wavefront, title_str = None):
        """ Plots an image on the global mask 
        given a wavefront """
        
        # full_img = np.reshape(wavefront, np.shape(self.global_mask))
        # img = np.ma.masked_array(full_img, self.global_mask)
        
        img = np.zeros(np.size(self.global_mask))
        flat_mask = self.global_mask.flatten()
        img[~flat_mask] = wavefront
        img = np.reshape(img, np.shape(self.global_mask))
        img = np.ma.masked_array(img, self.global_mask)

        plt.figure()
        plt.imshow(img, origin = 'lower', cmap='inferno')
        if title_str is not None:
            plt.title(title_str)
            
        plt.colorbar()
        
        
    def segment_scramble(self, mode_amp = 10e-6 ):
        """ Defines an initial segment scramble
        returning a masked array image """
        
        n_hex = n_hexagons(self.n_rings)
        # hex_data_len = np.sum(1-self.local_mask)
        
        file_path = self.savepath + 'segment_scramble.fits'
        try:
            masked_img = read_fits(file_path, is_ma = True)
            
            flat_img = masked_img.data.flatten()
            # Save values in segments' wavefront
            for k in range(n_hex):
                row_ids = self.valid_ids[k]
                self.segment[k].shape = flat_img[row_ids]
                
            self.shape = masked_img.data[~masked_img.mask]
        
        except FileNotFoundError:
            pass
        
        # Retrieve number of modes from the interaction matrix
        n_modes = int(np.shape(self.IM)[1]/n_hex)
        
        # Generate random mode coefficients
        mode_vec = np.random.randn(n_hex*n_modes)
        
        # Probability inversely proportional to spatial frequency
        m = int(np.ceil((np.sqrt(8*n_modes)-1.)/2.))
        freq_vec = np.repeat(np.arange(m)+1,np.arange(m)+1)
        prob_vec = 1./freq_vec[:n_modes]
        prob_vec_rep = np.tile(prob_vec,n_hex)
        
        # Modulate on the probability
        mode_vec = mode_vec * prob_vec_rep
        
        # Amplitude
        mode_vec *= mode_amp
        
        # Matrix product
        flat_img = self.IM*mode_vec
        
        # # Global modes
        # n_glob_modes = np.shape(self.glob_int_mat)[1]
        # glob_mode_vec = np.random.randn(n_glob_modes)
        # flat_img += self.glob_int_mat*glob_mode_vec
        
        # Save values in segments' wavefront
        for k in range(n_hex):
            row_ids = self.valid_ids[k]
            self.segment[k].shape = flat_img[row_ids]
        
        # Reshape and mask image
        img = np.zeros(np.size(self.global_mask))
        flat_mask = self.global_mask.flatten()
        img[~flat_mask] = flat_img
        img = np.reshape(img, np.shape(self.global_mask))
        masked_img = np.ma.masked_array(img, self.global_mask)
    
        # Save to fits
        write_to_fits(masked_img, file_path)
        self.shape = masked_img.data[~masked_img.mask]
        
        
    def compute_global_interaction_matrix(self, n_modes):
        """ Computes the global interaction matrix: 
            [n_pixels,n_modes] """
            
        int_mat_shape = [np.sum(1-self.global_mask),n_modes]
            
        file_path = self.savepath + str(n_modes) + 'modes_global_interaction_matrix.fits'
        try:
            self.glob_IM = read_fits(file_path, sparse_shape = int_mat_shape)
            return
        except FileNotFoundError:
            pass
        
        data_len = np.sum(1-self.global_mask)
        modes_data = np.zeros([data_len*n_modes])
        
        for j in range(n_modes):
            modes_data[data_len*j:data_len*(j+1)] = czern(j+1, self.global_mask)
            
        valid_ids = np.arange(data_len)
        row_indices = np.tile(valid_ids,n_modes)
        row = row_indices.flatten()
        
        mode_indices = np.arange(n_modes)
        col = np.repeat(mode_indices,data_len)
    
        self.glob_IM = csr_matrix((modes_data, (row,col)),  
                                  int_mat_shape)
        
        # Save to fits
        write_csr_to_fits(self.glob_IM, file_path)
        
        
        
    def compute_interaction_matrix(self, n_modes):
        
        hex_data_len = np.sum(1-self.local_mask)
        file_name = self.savepath + str(n_modes) + 'modes_interaction_matrix.fits'
        local_modes = np.zeros([hex_data_len*n_modes])
        for j in range(n_modes):
            local_modes[hex_data_len*j:hex_data_len*(j+1)] = czern(j+1, self.local_mask)
            
        self.IM = self._distribute_local_to_global(local_modes, file_name)
        
        # Distribute data to single segments
        n_hex = n_hexagons(self.n_rings)
        for k in range(n_hex):
            row_ids = self.valid_ids[k]
            self.segment[k].IM = self.IM[row_ids,n_modes*k:n_modes*(k+1)].toarray()
            
            
    def assemble_IFF_and_R_matrices(self):
        
        # Local matrices
        local_IFF_path = self.savepath + 'local_influence_functions_matrix.fits'
        local_R_path = self.savepath + 'local_reconstructor_matrix.fits'
        
        # Read/compute local IFF
        try:
            local_IFF = read_fits(local_IFF_path)
        except FileNotFoundError:
            self.LMC.simulate_influence_functions()
            local_IFF = self.LMC.sim_IFF
            write_to_fits(local_IFF, local_IFF_path)
        
        # Read/compute local Reconstructor
        try:
            Rec = read_fits(local_R_path)
        except FileNotFoundError:
            Rec = self.LMC.compute_reconstructor(local_IFF) 
            write_to_fits(Rec, local_R_path)
        
        for k in range(len(self.hex_centers)):
            self.segment[k].IFF = local_IFF
            self.segment[k].R = Rec
         
        # Assemble global matrices
        IFF_file_name = self.savepath + 'global_influence_functions_matrix.fits'
        R_file_name = self.savepath + 'global_reconstructor_matrix.fits'
        
        self.IFF = self._distribute_local_to_global(local_IFF.flatten(order='F'), IFF_file_name)
        self.R = self._distribute_local_to_global(Rec, R_file_name)
        
        
    def _distribute_local_to_global(self, local_data, file_path):
        
        hex_data_len = np.sum(1-self.local_mask)
        n_hex = n_hexagons(self.n_rings)
        
        data_shape = np.shape(local_data)
            
        N = int(np.size(local_data)/hex_data_len)
        glob_data_len = np.sum(1-self.global_mask)
        mat_shape = [glob_data_len,N*n_hex]
        
        if len(data_shape) > 1: #  local_data is a matrix
            if data_shape[0] < hex_data_len: # Reconstructor [Nacts,Npix]
                mat_shape = [N*n_hex, glob_data_len]
            # else: # IFF [Npix,Nacts]
            #     local_data = local_data.T
            local_data = local_data.flatten()
                  
        try:
            mat = read_fits(file_path, sparse_shape = mat_shape)
            return mat
        except FileNotFoundError:
            pass
        
        # # # row_indices = np.tile(self.hex_valid_ids, N)
        # # valid_ids = np.arange(glob_data_len)
        # # row_indices = np.tile(valid_ids, N)
        # valid_ids = np.arange(hex_data_len)
        # row_indices = np.tile(valid_ids, N*n_hex)
        row_indices = np.tile(self.valid_ids, N)
        row = row_indices.flatten()
        
        val_indices = np.arange(int(N*n_hex))
        col = np.repeat(val_indices, hex_data_len)
        
        data = np.tile(local_data, n_hex)
        
        if data_shape[0] < hex_data_len: # e.g. for the reconstructor [Nacts,Npix]
            mat = csr_matrix((data, (col,row)), mat_shape)
        else:
            mat = csr_matrix((data, (row,col)), mat_shape)
        
        # Save to fits
        write_csr_to_fits(mat, file_path)
        
        return mat
        
        
        
    def _assemble_global_mask(self):
        """ Assembles the global segmented mask from
        the local mask data """
        
        mask_file_path = self.savepath + 'global_mask.fits'
        ids_file_path = self.savepath + 'valid_ids.fits'
        try:
            self.global_mask = read_fits(mask_file_path, is_bool = True)
            self.valid_ids = read_fits(ids_file_path)
            return
        except FileNotFoundError:
            pass
        
        # Height of hex + gap
        L = self.gap + 2.*self.hex_side_len*SIN60
        
        # Full mask dimensions
        Ny = np.ceil((L*self.n_rings + self.hex_side_len*SIN60)*2*self.pix_scale)
        Nx = np.ceil((L*self.n_rings*SIN60 + self.hex_side_len*(0.5+COS60))*2*self.pix_scale)
        mask_shape = np.array([int(Ny),int(Nx)])
        mask = np.zeros(mask_shape, dtype = bool)
        
        # Hexagon centers pixel coordinates
        pix_coords = self.hex_centers*self.pix_scale + np.array([Nx,Ny])/2.
        
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
        
        hex_data_len = np.sum(1-self.local_mask)
        rep_pix_coords = np.repeat(pix_coords, hex_data_len, axis = 0)
        
        global_row_idx = (rep_local_row + rep_pix_coords[:,1]).astype(int)
        global_col_idx = (rep_local_col + rep_pix_coords[:,0]).astype(int)
        
        # Save valid hexagon indices
        hex_idx = global_col_idx + global_row_idx*int(Nx)
        hex_valid_ids = np.reshape(hex_idx,[n_hex,hex_data_len])

        # Data
        data = np.ones(len(rep_pix_coords), dtype=bool)
        
        # Sparse matrix assembly
        sparse_mat = csr_matrix((data, (global_row_idx, global_col_idx)),  
                                  mask_shape, dtype=bool)
        mask += sparse_mat
            
        self.global_mask = np.array(~mask)
        
        # Save valid hexagon indices
        global_data_len = np.sum(1-self.global_mask)
        flat_ids = np.arange(global_data_len)
        flat_img = np.zeros(np.size(self.global_mask))
        flat_mask = (self.global_mask.copy()).flatten()
        flat_img[~flat_mask] = flat_ids
        # flat_hex_ids = hex_valid_ids.flatten()
        self.valid_ids = (flat_img[hex_valid_ids]).astype(int)
        
        # Save to fits
        write_to_fits((self.global_mask).astype(np.uint8), mask_file_path)
        write_to_fits(self.valid_ids, ids_file_path)
        
        
        
    def _define_segment_array(self):
        """ Builds an array of Segment class objects,
        containing their center coordinates and the
        actuator positions """
        
        self.segment = []
        
        for k,coords in enumerate(self.hex_centers):
            self.segment.append(Segment(k, coords, self.LMC))
        
        
    def _initialize_global_actuator_coordinates(self):
        
        n_hex = n_hexagons(self.n_rings)
        local_coords = self.LMC.local_act_coords
        n_acts = len(self.LMC.local_act_coords)
        
        x_act_coords = np.repeat(local_coords[:,0], n_hex)
        x_hex_centers = np.tile(self.hex_centers[:,0], n_acts)
        x_act_coords += x_hex_centers
        y_act_coords = np.repeat(local_coords[:,1], n_hex)
        y_hex_centers = np.tile(self.hex_centers[:,1], n_acts)
        y_act_coords += y_hex_centers
        act_coords = np.vstack([x_act_coords, y_act_coords])
        self.act_coords = act_coords.T
        
        self.act_pos = np.zeros(n_acts*n_hex)
        
        
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
        
        
        