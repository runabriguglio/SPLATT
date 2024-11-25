import numpy as np
# import matplotlib.pyplot as plt

from deformable_mirror import DeformableMirror as DM

class Segment(DM):
    """ Class defining the single segment
    with center coordinates and actuator displacements """
    
    def __init__(self, center_coordinates, mask, valid_ids):
        
        self.center = center_coordinates
        self.mask = mask
        self.valid_ids = valid_ids
    
    
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
        
        if np.abs(len(new_act_coords)-len(self.act_coords)):
            raise NotImplementedError('Change in actuator number not implemented')
        
        # Actuator data
        self.act_coords = new_act_coords
        
        # Update IFF and R ...
        
    
    
    
        