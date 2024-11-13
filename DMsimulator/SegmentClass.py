import numpy as np
import matplotlib.pyplot as plt



class Segment():
    """ Class defining the single segment
    with center coordinates and actuator displacements """
    
    def __init__(self, segment_id, center_coordinates, LMC):
        
        self.id = segment_id
        self.center = center_coordinates
        
        self.LMC = LMC
        self.mask = LMC.local_mask
        self.shape = np.zeros(np.sum(1-self.mask))
        
        self.act_coords = LMC.local_act_coords
        self.act_pos = np.zeros(len(self.act_coords))
        
        
    def wavefront(self):
        """ Plots an image of the segment's current shape 
        on the local segment mask """
        
        # Project wavefront on mask
        mask = self.mask.copy()
        flat_mask = mask.flatten()
        image = np.zeros(np.size(flat_mask))
        image[~flat_mask] = self.shape
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
        """ Command the position of the segments's 
        actuators in absolute (default) or relative terms"""
        
        old_pos = self.act_pos
        
        if absolute_pos:
            new_pos = pos_cmd
        else:
            new_pos = pos_cmd + old_pos

        self.act_pos = new_pos
        self.shape += self.IFF @ (new_pos-old_pos)
        
        return new_pos
    
    
    def flatten(self):
        """ Computes and applies the flat command
        for the current segment shape """
        
        act_cmd = self._apply_shape(-self.shape)
        flat_rms = np.std(self.shape)
        
        return act_cmd, flat_rms
    
    
    def mirror_command(self, mode_amps):
        """ Computes and applies a modal amplitude
        command to the segment """
        
        # length check
        n_modes = np.shape(self.IM)[1]
        if len(mode_amps) < n_modes:
            amps = np.zeros(n_modes)
            amps[:len(mode_amps)] = mode_amps
            mode_amps = amps
        
        shape = self.IM @ mode_amps
        self._apply_shape(shape)
        
        fitting_err = np.std(self.shape - shape)
        
        return fitting_err
    
    def update_act_coords(self, new_act_coords):
        
        # Actuator data
        self.act_coords = new_act_coords
        self.act_pos = np.zeros(len(self.act_coords))
        
        self.IFF = self.LMC.simulate_influence_functions(new_act_coords)
        self.R = np.linalg.pinv(self.IFF)
        
        
    def _apply_shape(self, shape):
        """ Computes and applies a shape 
        command to the segment"""
        
        # Compute cmd and wavefront change
        act_cmd = self.R @ shape
        delta_wf = self.IFF @ act_cmd
        
        # Update act position and wavefront
        self.act_pos += act_cmd
        self.shape += delta_wf
        
        return act_cmd
    
    
    
        