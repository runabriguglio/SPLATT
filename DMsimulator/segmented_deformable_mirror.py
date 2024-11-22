import numpy as np
# import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
# from scipy.sparse import bsr_matrix

from mirror_segment import Segment
from hexagonal_geometry import HexGeometry
import matrix_calculator as matcalc
import read_and_write_fits as rwf

from deformable_mirror import DeformableMirror as DM

class SegmentedMirror(DM):
    """ Class defining the segmented deformable mirror """

    def __init__(self, TN):
        
        # Define geometric quantities
        self.geom = HexGeometry(TN)
        
        # Define mask and shape
        self.mask = self.geom.global_mask
        self.shape = np.zeros(np.sum(1-self.mask))
        
        # Initialize segmented array
        self._define_segment_array()
        self._initialize_actuator_coordinates()
        
        
    def compute_global_zern_matrix(self, n_modes):
        """ Computes or reads the globl Zernike modes interaction matrix """
        
        file_path = self.geom.savepath + str(n_modes) + 'modes_global_zernike_mat.fits'
        
        try:
            self.glob_ZM = rwf.read_fits(file_path)
            return
        except FileNotFoundError:
            self.glob_ZM = matcalc.compute_zernike_matrix(self.mask, n_modes)
            rwf.write_to_fits(self.glob_ZM, file_path)
            
        
    def compute_local_zern_matrix(self, n_modes):
        """ Computes or reads the local Zernike modes interaction matrix
        distributing the result to all segments and assemling the global
        sparse interaction matrix """
        
        file_path = self.geom.savepath + str(n_modes) + 'modes_local_zernike_mat.fits'
        sparse_file_path = self.geom.savepath + str(n_modes) + 'modes_global_sparse_zernike_mat.fits'
        
        try:
            local_ZM = rwf.read_fits(file_path)
        except FileNotFoundError:
            local_ZM = matcalc.compute_zernike_matrix(self.geom.local_mask, n_modes)
            rwf.write_to_fits(local_ZM, file_path)
        
        for k in range(self.geom.n_hex):
            self.segment[k].ZM = local_ZM
        
        self.ZM = self._distribute_local_to_global(local_ZM, sparse_file_path)
        
        
    def initialize_IFF_and_R_matrices(self, simulate:bool = True):
        """ Initializes the local IFF and R matrices (same for all segments)
        distributing the result to all segments and assemling the global
        sparse IFF and R matrices """
        
        # Influence functions
        local_IFF_path = self.geom.savepath + 'local_influence_functions_matrix.fits'
        global_IFF_path = self.geom.savepath + 'global_influence_functions_matrix.fits'
        
        # Read/compute local IFF
        try:
            local_IFF = rwf.read_fits(local_IFF_path)
        except FileNotFoundError:
            ref_act_coords = self.segment[0].act_coords
            if simulate:
                local_IFF = matcalc.simulate_influence_functions(ref_act_coords, self.geom.local_mask, self.geom.pix_scale)
            else:
                local_IFF = matcalc.calculate_influence_functions(ref_act_coords, self.geom.local_mask)
                rwf.write_to_fits(local_IFF, local_IFF_path)
        
        # Reconstructor
        local_R_path = self.geom.savepath + 'local_reconstructor_matrix.fits'
        global_R_path = self.geom.savepath + 'global_reconstructor_matrix.fits'
        
        # Read/compute local R
        try:
            local_R = rwf.read_fits(local_R_path)
        except FileNotFoundError:
            local_R = matcalc.compute_reconstructor(local_IFF)
            rwf.write_to_fits(local_R, local_R_path)
    
        # Distribute to all segments
        for k in range(self.geom.n_hex):
            self.segment[k].IFF = local_IFF
            self.segment[k].R = local_R
            
        self.IFF = self._distribute_local_to_global(local_IFF, global_IFF_path)
        self.R = self._distribute_local_to_global(local_R, global_R_path)
    
        
        
        
    def update_global_act_coords(self):
        """ Reads and saves the positions and (global)
        coordinates of all actuatrs"""
        
        # Initialize actuator positions and coordinates
        pos = self.segment[0].act_pos
        coords = self.segment[0].act_coords
        
        for k in range(1,self.geom.n_hex):
            pos = np.vstack([pos, self.segment[k].act_pos])
            glob_coords = self.segment[k].act_coords + self.segment[k].center
            coords = np.vstack([coords, glob_coords])
            
        self.act_pos = pos
        self.act_coords = coords
        
        
    def _define_segment_array(self):
        """ Builds an array of Segment class objects,
        containing their center coordinates and the
        actuator positions """
        
        self.segment = []
        
        for k,coords in enumerate(self.geom.hex_centers):
            self.segment.append(Segment(coords, self.geom.local_mask))
            
            
    def _initialize_actuator_coordinates(self):
        """ Initializes the local actuator coordinates 
        for all segments in the mirror """
        
        local_act_coords = self.geom.initialize_segment_act_coords()
        
        for k in range(self.geom.n_hex):
            self.segment[k].update_act_coords(local_act_coords)
            
        self.update_global_act_coords()

                
    def _distribute_local_to_global(self, local_data, file_path):
        """ Function to distribute the local_data matrix
        from the local mask to the global one """
        
        hex_data_len = np.sum(1-self.geom.local_mask)
        glob_data_len = np.sum(1-self.mask)
        
        N = int(np.size(local_data)/hex_data_len)
        
        mat_shape = [glob_data_len,N*self.geom.n_hex]
        
        data_shape = np.shape(local_data)
        
        if len(data_shape) > 1: #  local_data is a matrix
            if data_shape[0] < hex_data_len: # Reconstructor [Nacts,Npix]
                mat_shape = [N*self.geom.n_hex, glob_data_len]
                local_data = local_data.flatten()
            else: # IFF [Npix,Nacts]
                local_data = local_data.flatten(order='F')
                  
        try:
            mat = rwf.read_fits(file_path, sparse_shape = mat_shape)
            return mat
        except FileNotFoundError:
            pass
        
        row_indices = np.tile(self.geom.valid_ids, N)
        row = row_indices.flatten()
        
        val_indices = np.arange(int(N*self.geom.n_hex))
        col = np.repeat(val_indices, hex_data_len)
        
        data = np.tile(local_data, self.geom.n_hex)
        
        if data_shape[0] < hex_data_len: # e.g. for the reconstructor [Nacts,Npix]
            mat = csr_matrix((data, (col,row)), mat_shape)
        else:
            mat = csr_matrix((data, (row,col)), mat_shape)
        
        # Save to fits
        rwf.write_csr_to_fits(mat, file_path)
        
        return mat
        
    

        
            
    # def assemble_IFF_and_R_matrices(self, simulated_IFFs = False):
        
    #     # Local matrices
    #     local_IFF_path = self.savepath + 'local_influence_functions_matrix.fits'
    #     local_R_path = self.savepath + 'local_reconstructor_matrix.fits'
        

        
    #     # Read/compute local Reconstructor
    #     try:
    #         Rec = read_fits(local_R_path)
    #     except FileNotFoundError:
    #         Rec = self.LMC.compute_reconstructor(local_IFF) 
    #         write_to_fits(Rec, local_R_path)
        

    #         self.segment[k].R = Rec
         
    #     # Assemble global matrices
    #     IFF_file_name = self.savepath + 'global_influence_functions_matrix.fits'
    #     R_file_name = self.savepath + 'global_reconstructor_matrix.fits'
        
    #     self.IFF = self._distribute_local_to_global(local_IFF, IFF_file_name)
    #     self.R = self._distribute_local_to_global(Rec, R_file_name)
        
        
    # def _distribute_local_to_global(self, local_data, file_path):
    #     """ Function to distribute the local_data matrix
    #     from the local mask to the global one """
        
    #     hex_data_len = np.sum(1-self.local_mask)
    #     glob_data_len = np.sum(1-self.global_mask)
        
    #     n_hex = n_hexagons(self.n_rings)
    #     N = int(np.size(local_data)/hex_data_len)
        
    #     mat_shape = [glob_data_len,N*n_hex]
        
    #     data_shape = np.shape(local_data)
    #     if len(data_shape) > 1: #  local_data is a matrix
    #         if data_shape[0] < hex_data_len: # Reconstructor [Nacts,Npix]
    #             mat_shape = [N*n_hex, glob_data_len]
    #             local_data = local_data.flatten()
    #         else: # IFF [Npix,Nacts]
    #             local_data = local_data.flatten(order='F')
                  
    #     try:
    #         mat = rwf.read_fits(file_path, sparse_shape = mat_shape)
    #         return mat
    #     except FileNotFoundError:
    #         pass
                
    #     row_indices = np.tile(self.valid_ids, N)
    #     row = row_indices.flatten()
        
    #     val_indices = np.arange(int(N*n_hex))
    #     col = np.repeat(val_indices, hex_data_len)
        
    #     data = np.tile(local_data, n_hex)
        
    #     if data_shape[0] < hex_data_len: # e.g. for the reconstructor [Nacts,Npix]
    #         mat = csr_matrix((data, (col,row)), mat_shape)
    #     else:
    #         mat = csr_matrix((data, (row,col)), mat_shape)
        
    #     # Save to fits
    #     rwf.write_csr_to_fits(mat, file_path)
        
    #     return mat
        
        
        