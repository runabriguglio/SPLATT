import numpy as np
import matplotlib.pyplot as plt

from matrix_calculator import matmul

class DeformableMirror():
    
    def __init__(self):
        pass
    
            
    def surface(self, surf2plot = None, plot_title:str = None):
        """ Plots surf2plot or (default) the segment's
        current shape on the DM mask """
        
        if surf2plot is None:
            surf2plot = self.shape
        
        # Project wavefront on mask
        mask = self.mask.copy()
        flat_mask = mask.flatten()
        image = np.zeros(np.size(flat_mask))
        image[~flat_mask] = surf2plot
        image = np.reshape(image, np.shape(mask))
        image = np.ma.masked_array(image, mask)
        
        # Plot image
        plt.figure()
        plt.imshow(image, origin = 'lower', cmap = 'inferno')
        plt.colorbar()
        
        if plot_title is not None:
            plt.title(plot_title)
        
        
    def get_position(self):
        """ Wrapper to read the current 
        position of the segments's actuators """
        
        pos = self.act_pos
        x,y = self.act_coords[:,0], self.act_coords[:,1] 
        
        plt.figure()
        plt.scatter(x,y, c=pos, s=100, cmap='inferno')
        plt.colorbar()
        
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
            
        act_cmd = matmul(self.R,delta_shape)
        self.mirror_command(act_cmd)
    
    
    def mirror_command(self, cmd_amps, zonal:bool = True):
        """ Computes and applies a zonal/modal amplitude
        command to the segment """
        
        if zonal is False: # convert modal to zonal
            
            mode_amps = cmd_amps.copy() # rename
            
            # length check
            n_modes = np.shape(self.ZM)[1]
            if len(mode_amps) < n_modes:
                amps = np.zeros(n_modes)
                amps[:len(mode_amps)] = mode_amps
                mode_amps = amps
        
            shape = matmul(self.ZM,mode_amps)
            cmd_amps = matmul(self.R,shape)
            
        # Update actuator position
        self.act_pos += cmd_amps
        
        # Update mirror shape
        delta_shape = matmul(self.IFF,cmd_amps)
        self.shape += delta_shape

