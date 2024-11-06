import numpy as np
import os

from read_config import readConfig
from rotate_coordinates import cw_rotate as crot
from read_and_write_fits import write_to_fits
from read_and_write_fits import read_fits 

# Useful variables
SIN60 = np.sin(np.pi/3.)
COS60 = np.cos(np.pi/3.)

def n_hexagons(n_rings):
    return int(1 + (6 + n_rings*6)*n_rings/2)


class Segment():
    """ Class defining the single segment
    with center coordinates and actuator displacements """
    
    def __init__(self, center_coordinates):
        
        self.center = center_coordinates
        self.act_pos = []
        

class DM():
    """ Class defining the segmented deformable mirror """

    def __init__(self, TN):
        
        # Read configuration files
        config_path = './config_par_' + TN + '.yaml'
        dm_par, opt_par = readConfig(config_path)
        
        self.gap = dm_par[0]
        self.hex_side_len = dm_par[1]
        self.n_rings = int(dm_par[2])

        self.pix_scale = opt_par[0]
        
        self.savepath = './' + TN + '/'
        self.n_hex = n_hexagons(self.n_rings)
        
        try: # Create new folder
            os.mkdir(TN)
        except FileExistsError: # Folder already existing
            pass
        
        self._compute_segment_centers()
        self._define_segment_array()
        
        
        
    def _define_segment_array(self):
        """ Builds an array of Segment class objects,
        containing their center coordinates and the
        actuator positions """
        
        self.segments = []
        
        for k,coords in enumerate(self.hex_centers):
            self.segments.append(Segment(coords, n))
        
        
    def _compute_segment_centers(self):
        """ Defines and saves the coordinates of the 
        centers of all hexagonal segments """
        
        file_path = self.savepath + 'hex_centers_coords.fits'
        try:
            self.hex_centers = read_fits(file_path)
            return
        except FileNotFoundError:
            pass
        
        # Number of hexes
        self.hex_centers = np.zeros([n_hexagons(self.n_rings),2])
        
        # Height of hex + gap
        L = self.gap + 2.*self.hex_side_len*SIN60
        angles = np.pi/3. * (np.arange(5)+1)
        
        for ring_ctr in range(self.n_rings):
            
            hex_ctr = n_hexagons(ring_ctr)
            ring_ctr += 1
            R = L*ring_ctr
            
            aux = np.array([R*SIN60,R*COS60])
            self.hex_centers[hex_ctr,:] = aux
            
            self.hex_centers[hex_ctr+ring_ctr:hex_ctr+6*ring_ctr:ring_ctr,:] = crot(aux, angles)
            
            if ring_ctr > 1:
                for j in range(ring_ctr-1):
                    shift = self.gap + 2.*self.hex_side_len*SIN60
                    aux[0] = self.hex_centers[hex_ctr,0] 
                    aux[1] = self.hex_centers[hex_ctr,1] - (j+1)*shift
                    self.hex_centers[hex_ctr+j+1,:] = aux
                    self.hex_centers[hex_ctr+j+1+ring_ctr:hex_ctr+j+1+6*ring_ctr:ring_ctr,:] = crot(aux, angles)

        # Save to fits
        write_to_fits(self.hex_centers, file_path)
        