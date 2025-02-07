import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

from matrix_calculator import matmul

class DeformableMirror():
    """ 
    This class defines the methods for a generic deformable mirror.
    It implements the following functions:
        - surface(): plot the current shape or a shape in input
        - get_position(): plot the actuator position at the actuaor coordinates
        - apply_flat(): apply a flat command of the current shape with
                        an (optional) offset
        - mirror_command(): apply a zonal or modal command to the DM
    """
    
    def __init__(self):
        pass

    def surface(self, surf2plot = None, plt_title:str = None):
        """
        Plots surf2plot or (default) the segment's
        current shape on the DM mask

        Parameters
        ----------
        surf2plot : ndarray(float) [Npix], optional
            The shape to plot. The default (no input given)
            is the DM current shape.
        plt_title : str, optional
            The plot title. The default is None.

        Returns
        -------
        None.

        """
        
        if surf2plot is None:
            surf2plot = self.shape
        
        # Project wavefront on mask
        mask = self.mask.copy()
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
            plt_title = f"RMS: {np.std(self.shape):.2e}"
            
        # plt.axis([150,200,150,250]) # debug
        plt.axis('off')
        plt.title(plt_title)

    def get_position(self, act_pix_size:float = 10):
        """
        Readsa and plots the current actuators' 
        position at the actuators' coordinates 

        Parameters
        ----------
        act_pix_size : float, optional
            The actuator size in pixels. The default is 10.

        Returns
        -------
        pos : ndarray(float) [Nacts,]
            Array of current actuator positions on the segment.

        """
        
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
        offset : ndarray(float) [Npix], optional
            The shape offset from which to compute the flat.
            The default is no offset.

        Returns
        -------
        act_cmd : ndarray(float) [Nacts]
            The actoator position command to achieve the flat.
        flat_rms : float
            The standard deviation of the obtained shape.

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
        """
        Computes and applies a zonal/modal amplitude
        command to the segment

        Parameters
        ----------
        cmd_amps : ndarray(float) [Nacts]
            The array of zonal (modal) amplitudes.
            
        absolute_delta_pos : bool, optional
            True if the command is absolute and not
            relative to the current actuators' position.
            The default is False.
            
        modal : bool, optional
            True for modal amplitudes. The default is False.

        Returns
        -------
        None.

        """
        
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

