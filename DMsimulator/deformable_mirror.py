import numpy as np
import matplotlib.pyplot as plt

from matrix_calculator import matmul

class DeformableMirror():
    """ This class defines the methods for a generic deformable mirror"""
    
    def __init__(self):
        pass
    
            
    def surface(self, surf2plot = None, plt_title:str = None):
        """ Plots surf2plot or (default) the segment's
        current shape on the DM mask """
        
        if surf2plot is None:
            surf2plot = self.shape
        
        # Project wavefront on mask
        mask = self.mask.copy()
        
        # if is_global:
        #     pix_ids = ~mask.flatten()
        # else:
        pix_ids = self.valid_ids.flatten()
            
        image = np.zeros(np.size(mask))
        image[pix_ids] = surf2plot
        image = np.reshape(image, np.shape(mask))
        image = np.ma.masked_array(image, mask)
        
        # Plot image
        plt.figure()
        plt.imshow(image, origin = 'lower', cmap = 'inferno')
        plt.colorbar()
        
        if plt_title is None:
            plt_title = 'RMS: ' + str(np.std(self.shape))
            
        plt.title(plt_title)
        
        
    def get_position(self, act_pix_size:float = 10):
        """ Wrapper to read the current 
        position of the segments's actuators """
        
        pos = self.act_pos.copy()
        x,y = self.act_coords[:,0], self.act_coords[:,1] 
        
        plt.figure()
        plt.scatter(x,y, c=pos, s=act_pix_size**2, cmap='inferno')
        plt.colorbar()
        
        return pos
    
    
    def apply_flat(self, offset = None):
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
            
        act_cmd = matmul(self.R, delta_shape)
        self.mirror_command(act_cmd)
        
        flat_rms = np.std(self.shape)
        
        if offset is not None:
            flat_rms = np.std(self.shape - offset)
        
        return act_cmd, flat_rms
    
    
    def mirror_command(self, cmd_amps, absolute_delta_pos:bool = False, modal:bool = False):
        """ Computes and applies a zonal/modal amplitude
        command to the segment """
        
        if modal: # convert modal to zonal
            
            mode_amps = cmd_amps.copy() # rename
            
            # length check
            n_modes = np.shape(self.ZM)[1]
            if len(mode_amps) < n_modes:
                amps = np.zeros(n_modes)
                amps[:len(mode_amps)] = mode_amps
                mode_amps = amps
        
            shape = matmul(self.ZM, mode_amps)
            cmd_amps = matmul(self.R, shape)
        
        # Position (zonal) command
        delta_pos = cmd_amps - absolute_delta_pos * self.act_pos

        # Update positions and shape
        self.act_pos += delta_pos
        self.shape += matmul(self.IFF, delta_pos)

