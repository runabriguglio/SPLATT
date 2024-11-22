import numpy as np
# import matplotlib.pyplot as plt

from deformable_mirror import DeformableMirror as DM

class Segment(DM):
    """ Class defining the single segment
    with center coordinates and actuator displacements """
    
    def __init__(self, center_coordinates, mask):
        
        self.center = center_coordinates
        
        self.mask = mask
        self.shape = np.zeros(np.sum(1-self.mask))
    
    
    def update_act_coords(self, new_act_coords):
        """
        Updates the actuator coordinates to new_act_coordinates 
        (resetting their position to zero) 
        

        Parameters
        ----------
        new_act_coords : ndarray([n_acts,2])
            Vector containing the new actuator coordinates.

        Returns
        -------
        None.

        """
        
        # Actuator data
        self.act_coords = new_act_coords
        self.act_pos = np.zeros(len(self.act_coords))
        
    
    
    
        