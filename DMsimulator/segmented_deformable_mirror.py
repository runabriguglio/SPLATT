import numpy as np
# import matplotlib.pyplot as plt
# from scipy.sparse import csr_matrix
# from scipy.sparse import bsr_matrix

from mirror_segment import Segment
# from hexagonal_geometry import HexGeometry
import matrix_calculator as matcalc
import read_and_write_fits as rwf

from deformable_mirror import DeformableMirror as DM

class SegmentedMirror(DM):
    """ Class defining the segmented deformable mirror """

    def __init__(self, mirror_geometry):
        
        # Define geometric quantities
        self.geom = mirror_geometry
        
        # Define mask and shape
        self.mask = self.geom.global_mask
        
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
            mask = self.mask.copy()
            scrambled_ZM = matcalc.compute_zernike_matrix(mask, n_modes)
            ids = np.zeros(np.size(mask),dtype=int)
            ids[~mask.flatten()] = np.arange(np.sum(1-mask))
            scrambled_ids = ids[self.valid_ids.flatten()]
            self.glob_ZM = scrambled_ZM[scrambled_ids,:]
            rwf.write_to_fits(self.glob_ZM, file_path)
            
        
    def compute_local_zern_matrix(self, n_modes):
        """ Computes or reads the local Zernike modes interaction matrix
        distributing the result to all segments and assemling the global
        sparse interaction matrix """
        
        file_path = self.geom.savepath + str(n_modes) + 'modes_local_zernike_mat.fits'
        # n_pix = np.sum(1-self.geom.local_mask)
        # n_hex = self.geom.n_hex
        # mat_shape = [n_pix*n_hex, n_modes*n_hex]
        
        try:
            ZM = rwf.read_fits(file_path)#, sparse_shape = mat_shape)
        except FileNotFoundError:
            loc_ZM = matcalc.compute_zernike_matrix(self.geom.local_mask, n_modes)
            ZM = np.tile(loc_ZM,(self.geom.n_hex,1,1))#self._distribute_local_to_global(loc_ZM)
            rwf.write_to_fits(ZM, file_path)#rwf.write_csr_to_fits(ZM, file_path)
        
        self.ZM = ZM
        for k,segment in enumerate(self.segment):
            segment.ZM = self.ZM[k]#(self.ZM[n_pix*k:n_pix*(k+1),n_modes*k:n_modes*(k+1)]).toarray()
        
        
    def initialize_IFF_and_R_matrices(self, simulate:bool = True):
        """ Initializes the local IFF and R matrices (same for all segments)
        distributing the result to all segments and assemling the global
        sparse IFF and R matrices """
        
        IFF_path = self.geom.savepath + 'global_influence_functions_matrix.fits'
        R_path = self.geom.savepath + 'global_reconstructor_matrix.fits'
        
        # n_pix = np.sum(1-self.geom.local_mask)
        # n_hex = self.geom.n_hex
        # n_acts = int(len(self.act_coords)/n_hex)
        
        # IFF_shape = [n_pix*n_hex, n_hex*n_acts]
        # R_shape = [n_acts*n_hex, n_hex*n_pix]
        
        # Read/compute IFF and R
        try:
            self.IFF = rwf.read_fits(IFF_path)#, sparse_shape = IFF_shape)
            self.R = rwf.read_fits(R_path)#, sparse_shape = R_shape)
            
        except FileNotFoundError:
            IFF_cube = self._compute_local_IFF_cube()
            
            loc_IFF = matcalc.cube2mat(IFF_cube)
            self.IFF = np.tile(loc_IFF,(self.geom.n_hex,1,1))#self._distribute_local_to_global(loc_IFF)
            rwf.write_to_fits(self.IFF, IFF_path)#rwf.write_csr_to_fits(self.IFF, IFF_path)
            
            loc_R = matcalc.compute_reconstructor(loc_IFF)
            self.R = np.tile(loc_R,(self.geom.n_hex,1,1))#self._distribute_local_to_global(loc_R)
            rwf.write_to_fits(self.R, R_path)#rwf.write_csr_to_fits(self.R, R_path)
    
        # Distribute to all segments
        for k,segment in enumerate(self.segment):
            segment.IFF = self.IFF[k]#(self.IFF[n_pix*k:n_pix*(k+1), n_acts*k:n_acts*(k+1)]).toarray()
            segment.R = self.R[k]#(self.R[n_acts*k:n_acts*(k+1), n_pix*k:n_pix*(k+1)]).toarray()
            
            
    def _compute_local_IFF_cube(self, segment_id:int = None, simulate:bool = True):
        """ Reads or computes the IFF image cube for the actuator
        coordinates of the segment with given segment_id"""
        
        if segment_id is None:
            file_name = 'REF'
            segment_id = 0
        else:
            file_name = 'segment' + str(segment_id)
        
        file_path = self.geom.savepath + file_name + '_IFF_image_cube.fits'
        
        if file_name == 'REF':
            try:
                IFF_cube = rwf.read_fits(file_path, is_ma = True)
                return IFF_cube
            except FileNotFoundError:
                pass
        
        ref_act_coords = self.segment[segment_id].act_coords
        if simulate:
            IFF_cube = matcalc.simulate_influence_functions(ref_act_coords, self.geom.local_mask, self.geom.pix_scale)
        else:
            IFF_cube = matcalc.calculate_influence_functions(ref_act_coords, self.geom.local_mask, self.geom.act_radius/self.geom.hex_side_len)
        rwf.write_to_fits(IFF_cube, file_path)
        
        return IFF_cube
            
            
    def update_act_coords(self, segment_id:int, new_act_coords):
        """ Updates the actuator coordinates of the segment with given
        segment_id, updating the IFF and R matrices accordingly """
        
        n_acts = int(len(self.act_coords)/self.geom.n_hex)
        n_new = len(new_act_coords) 
        
        if n_new > n_acts:
            raise NotImplementedError('Current number of actuators is ' + str(n_acts) + 
                                      ' an increase to ' + str(n_new) + ' is not yet implemented')
        elif n_new < n_acts:
            nan_padding = np.empty([n_acts - n_new,2])
            nan_padding.fill(np.nan)
            new_act_coords = np.vstack([new_act_coords, nan_padding])
            
        act_ids = np.arange(segment_id*n_acts,(segment_id+1)*n_acts)
        
        # Update actuator coordinates
        self.segment[segment_id].act_coords = new_act_coords
        self.act_coords[act_ids] = new_act_coords + self.segment[segment_id].center
        
        # Compute IFF and R matrices
        IFF_cube = self._compute_local_IFF_cube(segment_id, simulate=True)
        loc_IFF = matcalc.cube2mat(IFF_cube)
        loc_R = matcalc.compute_reconstructor(loc_IFF)
        
        # n_pix = np.sum(1-self.segment[segment_id].mask)
        # pix_ids = np.arange(segment_id*n_pix,(segment_id+1)*n_pix)
        
        # Update IFF and R matrices
        # using += in order to update global at the same time
        old_IFF = self.segment[segment_id].IFF
        self.segment[segment_id].IFF += loc_IFF - old_IFF 
        old_R = self.segment[segment_id].R
        self.segment[segment_id].R += loc_R - old_R
        # self.IFF[pix_ids, act_ids] = loc_IFF
        # self.R[act_ids, pix_ids] = loc_R
        
        # Reset shape and position
        self.segment[segment_id].shape *= 0
        self.segment[segment_id].act_pos *= 0
            
                
        
        
    def _define_segment_array(self):
        """ Builds an array of Segment class objects,
        containing their center coordinates and the
        actuator positions """
        
        loc_mask = self.geom.local_mask
        
        local_ids = np.arange(np.size(loc_mask))
        local_valid_ids = local_ids[~loc_mask.flatten()]
        
        self.valid_ids = self.geom.valid_ids
        
        # Initialize segment array
        self.segment = []
        
        for k,coords in enumerate(self.geom.hex_centers):
            self.segment.append(Segment(coords, self.geom.local_mask, local_valid_ids))
            
            
    def _initialize_actuator_coordinates(self):
        """ Initializes the local actuator coordinates 
        for all segments in the mirror """
        
        local_act_coords = self.geom.initialize_segment_act_coords()
        
        n_acts = len(local_act_coords)
        n_hex = self.geom.n_hex
        n_pix = np.sum(1-self.geom.local_mask)
        
        self.shape = np.zeros(np.sum(1-self.mask))
        self.act_coords = np.tile(local_act_coords,(n_hex,1))
        self.act_pos = np.zeros(n_hex*n_acts)
        
        for k, segment in enumerate(self.segment):
            segment.act_coords = local_act_coords
            segment.act_pos = self.act_pos[n_acts*k:n_acts*(k+1)]
            segment.shape = self.shape[n_pix*k:n_pix*(k+1)]
            self.act_coords[n_acts*k:n_acts*(k+1),:] += segment.center
            
            
    # def _distribute_local_to_global(self, local_mat):
    #     """ Function to distribute the local_data matrix
    #     from the local mask to the global one """
        
    #     mat_shape = np.shape(local_mat)
    #     n_hex = self.geom.n_hex
    #     sparse_shape = np.dot(n_hex,mat_shape)
        
    #     row, col = matcalc.get_sparse_ids(mat_shape, n_hex)
    #     local_data = local_mat.flatten(order = 'F')
            
    #     data = np.tile(local_data,n_hex)
        
    #     sparse_mat = csr_matrix((data, (row,col)), sparse_shape)
        
    #     return sparse_mat

        
        