import numpy as np
# import matplotlib.pyplot as plt

from mirror_segment import Segment
# from hexagonal_geometry import HexGeometry
import matrix_calculator as matcalc
# import read_and_write_fits as myfits
import my_fits_package as myfits

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
            self.glob_ZM = myfits.read_fits(file_path)
            return
        except FileNotFoundError:
            mask = self.mask.copy()
            scrambled_ZM = matcalc.compute_zernike_matrix(mask, n_modes)
            ids = np.zeros(np.size(mask),dtype=int)
            ids[~mask.flatten()] = np.arange(np.sum(1-mask))
            scrambled_ids = ids[self.valid_ids.flatten()]
            self.glob_ZM = scrambled_ZM[scrambled_ids,:]
            myfits.write_to_fits(self.glob_ZM, file_path)
            
        
    def compute_local_zern_matrix(self, n_modes):
        """ Computes or reads the local Zernike modes interaction matrix
        distributing the result to all segments and assemling the global
        sparse interaction matrix """
        
        file_path = self.geom.savepath + str(n_modes) + 'modes_local_zernike_mat.fits'
        
        try:
            ZM = myfits.read_fits(file_path)
        except FileNotFoundError:
            loc_ZM = matcalc.compute_zernike_matrix(self.geom.local_mask, n_modes)
            ZM = np.tile(loc_ZM,(self.geom.n_hex,1,1))
            myfits.write_to_fits(ZM, file_path)
        
        self.ZM = ZM
        for k,segment in enumerate(self.segment):
            segment.ZM = self.ZM[k]
            
        
    def initialize_IFF_and_R_matrices(self, simulate:bool = True):
        """ Initializes the local IFF and R matrices (same for all segments)
        distributing the result to all segments and assemling the global
        sparse IFF and R matrices """
        
        self.IFF_path = self.geom.savepath + 'global_influence_functions_matrix.fits'
        self.R_path = self.geom.savepath + 'global_reconstructor_matrix.fits'
        
        # Read/compute IFF and R
        try:
            self.IFF = myfits.read_fits(self.IFF_path)
            self.R = myfits.read_fits(self.R_path)
            
        except FileNotFoundError:
            IFF_cube = self._compute_local_IFF_cube(self.geom.local_act_coords)
            
            loc_IFF = matcalc.cube2mat(IFF_cube)
            self.IFF = np.tile(loc_IFF,(self.geom.n_hex,1,1))
            myfits.write_to_fits(self.IFF, self.IFF_path)
            
            loc_R = matcalc.compute_reconstructor(loc_IFF)
            self.R = np.tile(loc_R,(self.geom.n_hex,1,1))
            myfits.write_to_fits(self.R, self.R_path)
    
        # Distribute to all segments
        for k,segment in enumerate(self.segment):
            segment.IFF = self.IFF[k]
            segment.R = self.R[k]
            
            
    def _compute_local_IFF_cube(self, ref_act_coords, segment_id = None, simulate:bool = True):
        """ Reads or computes the IFF image cube for the actuator
        coordinates of the segment with given segment_id"""
        
        if segment_id is None:
            file_path = self.geom.savepath + 'Reference_IFF_image_cube.fits'
            try:
                IFF_cube = myfits.read_fits(file_path, is_ma = True)
                return IFF_cube
            except FileNotFoundError:
                pass   
        
        if simulate:
            IFF_cube = matcalc.simulate_influence_functions(ref_act_coords, self.geom.local_mask, self.geom.pix_scale)
        else:
            IFF_cube = matcalc.calculate_influence_functions(ref_act_coords, self.geom.local_mask, self.geom.act_radius/self.geom.hex_side_len)
       
        for k,idx in enumerate(segment_id):
            file_path = self.geom.savepath + 'segment' + str(idx) + '_IFF_image_cube.fits'
            myfits.write_to_fits(IFF_cube, file_path)
        
        return IFF_cube
            
            
    def update_act_coords(self, segment_ids, new_act_coords, do_save:bool = True):
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
            
        # Compute IFF and R matrices
        IFF_cube = self._compute_local_IFF_cube(new_act_coords, segment_ids, simulate=True)
        loc_IFF = matcalc.cube2mat(IFF_cube)
        loc_R = matcalc.compute_reconstructor(loc_IFF)
            
        for k,segm_id in enumerate(segment_ids):
            act_ids = np.arange(segm_id*n_acts,(segm_id+1)*n_acts)
            seg = self.segment[segm_id]
            
            # Update actuator coordinates
            old_act_coords = seg.act_coords# - self.segment[segment_id].center
            seg.act_coords += new_act_coords - old_act_coords
            self.act_coords[act_ids] = new_act_coords + seg.center
            
            # Update IFF and R matrices
            # using += in order to update global at the same time
            old_IFF = seg.IFF
            seg.IFF += loc_IFF - old_IFF 
            old_R = seg.R
            seg.R += loc_R - old_R
            
            # Reset shape and position
            seg.shape *= 0
            seg.act_pos *= 0
        
        # Save to .fits
        if do_save:
            myfits.write_to_fits(self.IFF, self.IFF_path)
            myfits.write_to_fits(self.R, self.R_path)
            myfits.write_to_fits(self.act_coords, self.coords_path)
            
                
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
        
        self.geom.local_act_coords = self.geom.initialize_segment_act_coords()
        n_acts = len(self.geom.local_act_coords)
        n_hex = self.geom.n_hex
        n_pix = np.sum(1-self.geom.local_mask)
        
        self.coords_path = self.geom.savepath + 'global_actuator_coordinates.fits'
        try:
            self.act_coords = myfits.read_fits(self.coords_path)
        except FileNotFoundError:
            self.act_coords = np.tile(self.geom.local_act_coords,(n_hex,1))
            self.act_coords += np.tile(self.geom.hex_centers,*1,n_acts()).reshape([n_hex*n_acts,2])
        
        self.shape = np.zeros(np.sum(1-self.mask))
        self.act_pos = np.zeros(len(self.act_coords))
        
        for k, segment in enumerate(self.segment):
            segment.act_coords = self.act_coords[n_acts*k:n_acts*(k+1),:] - segment.center
            segment.act_pos = self.act_pos[n_acts*k:n_acts*(k+1)]
            segment.shape = self.shape[n_pix*k:n_pix*(k+1)]
            
        # Save cooeds to fits
        myfits.write_to_fits(self.act_coords, self.coords_path)

        
        