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
        self.act_pix_size = LMC.act_radius*LMC.pix_scale
        
        
    def surface(self):
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
        x,y = self.act_coords[:,0], self.act_coords[:,1] 
        
        plt.figure()
        plt.scatter(x,y,c=pos, s=100, cmap='inferno')
        plt.colorbar()
        plt.title('Segment ' + str(self.id) + ' actuator command')
        
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
    
    
    def flatten(self, offset = None):
        """
        Computes and applies the flat command
        for the current segment shape minus a given offset

        Parameters
        ----------
        offset : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        act_cmd : TYPE
            DESCRIPTION.
        flat_rms : TYPE
            DESCRIPTION.

        """
        
        delta_shape = -self.shape
        
        if offset is not None:
            delta_shape += offset
            
        act_cmd = self.R @ delta_shape
        self.mirror_command(act_cmd)
    
    
    def mirror_command(self, cmd_amps, zonal:bool = True):
        """ Computes and applies a zonal/modal amplitude
        command to the segment """
        
        if zonal is False: # convert modal to zonal
            
            mode_amps = cmd_amps.copy() # rename
            
            # length check
            n_modes = np.shape(self.IM)[1]
            if len(mode_amps) < n_modes:
                amps = np.zeros(n_modes)
                amps[:len(mode_amps)] = mode_amps
                mode_amps = amps
        
            shape = self.IM @ mode_amps
            cmd_amps = self.R @ shape
            
        # Update actuator position
        self.act_pos += cmd_amps
        
        # Update mirror shape
        delta_shape = self.IFF @ cmd_amps
        self.shape += delta_shape
    
    
    def update_act_coords(self, new_act_coords, simulate:bool = True):
        """
        Updates the actuator coordinates to new_act_coordinates 
        (resetting their position to zero) and computes the 
        

        Parameters
        ----------
        new_act_coords : ndarray([n_acts,2])
            Vector containing the new actuator coordinates.
        simulate : bool, optional
            Variable for the simulation of the IFF
            rather than their (computationally expensive)
            FE simulation. The default is True.

        Returns
        -------
        None.

        """
        
        # Actuator data
        self.act_coords = new_act_coords
        self.act_pos = np.zeros(len(self.act_coords))
        
        if simulate:
            self.IFF = self.LMC.simulate_influence_functions(new_act_coords)
        else:
            self.IFF = self.LMC.compute_influence_functions(new_act_coords)
            
        self.R = np.linalg.pinv(self.IFF)
        
    
    
    
        