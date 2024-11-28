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
        
    
    
    
        