import numpy as np
import matplotlib.pyplot as plt

from read_config import readConfig
from zernike_polynomials import computeZernike as czern

from scipy.sparse import csr_matrix 


# Useful variables
SIN60 = np.sin(np.pi/3.)
COS60 = np.cos(np.pi/3.)


def cw_rotate(vec, angle):
    c = np.cos(angle)
    s = np.sin(angle)
    rot_mat = [[c,s],[-s,c]]
    
    if np.shape(vec)[0] > 2:
        aux_vec = rot_mat @ vec.transpose()
        rot_vec = aux_vec.transpose()
    else:
        rot_vec = rot_mat @ vec
    
    return rot_vec


def rotate_by_60deg(vec):
    # Wrapper to cw_rotate()
    cw_angle = np.pi/3.
    rot_vec = cw_rotate(vec, cw_angle)
    
    return rot_vec

def n_hexagons(n_rings):
    return int(1 + (6 + n_rings*6)*n_rings/2)


class Hexagons():

    def __init__(self, configfile):

        # Read configuration files
        path = './' + configfile
        par = readConfig(path)

        self.g = par[0]
        self.hex_side_len = par[1]
        self.n_rings = int(par[2])
        self.pix_scale = par[3]
        self.angle = par[4]
    

    def find_hex_coordinates(self):
        """ Defines and saves the coordinates of the 
        centers of all hexagonal segments """
        
        # Number of hexes
        self.hex_centers = np.zeros([n_hexagons(self.n_rings),2])
        
        # Height of hex + gap
        L = self.g + 2.*self.hex_side_len*SIN60
        
        for ring_ctr in range(self.n_rings):
            
            hex_ctr = n_hexagons(ring_ctr)
            ring_ctr += 1
            R = L*ring_ctr
            
            aux = np.array([R*SIN60,R*COS60])
            self.hex_centers[hex_ctr,:] = aux
            
            for i in range(5):
                aux = rotate_by_60deg(aux)
                self.hex_centers[hex_ctr+(i+1)*ring_ctr,:] = aux
            
            if ring_ctr > 1:
                for j in range(ring_ctr-1):
                    shift = self.g + 2.*self.hex_side_len*SIN60
                    aux[0] = self.hex_centers[hex_ctr,0] 
                    aux[1] = self.hex_centers[hex_ctr,1] - (j+1)*shift
                    self.hex_centers[hex_ctr+j+1,:] = aux
                    for i in range(5):
                        aux = rotate_by_60deg(aux)
                        self.hex_centers[hex_ctr+j+1+(i+1)*ring_ctr,:] = aux
                    
        # Print to a .fits file
        
        # return self.hex_centers


    def define_hex_mask(self):
        """ Forms the hexagonal mask to be placed 
        at the given hexagonal coordinates """
        L = self.hex_side_len*self.pix_scale
        
        # Consider points in the upper-left of the hexagon
        max_y = int(L*SIN60)
        max_x = int(L/2.+L*COS60)
        
        mask_ul = np.fromfunction(lambda i,j: j > L/2. + i/SIN60*COS60, [max_y,max_x])
        
        mask = np.zeros([2*max_y,2*max_x], dtype=bool)
        
        mask[0:max_y,max_x:] = mask_ul # upper left
        mask[0:max_y,0:max_x] = np.flip(mask_ul,1) # upper right
        mask[max_y:,:] = np.flip(mask[0:max_y,:],0) # lower
                    
        self.hex_mask = mask
        
        x = np.arange(2*max_x)
        y = np.arange(2*max_y)
        X,Y = np.meshgrid(x,y)
        valid_X = X[~self.hex_mask]
        valid_Y = Y[~self.hex_mask]
        
        self.valid_row = valid_Y - max_y
        self.valid_col = valid_X - max_x
        
        
    def define_segmented_mask(self):
        """ Forms the segment mask """
        
        # Height of hex + gap
        L = self.g + 2.*self.hex_side_len*SIN60
        
        # Full mask dimensions
        Ny = (L*self.n_rings + self.hex_side_len*SIN60)*2*self.pix_scale
        Nx = (L*self.n_rings*SIN60 + self.hex_side_len*(0.5+COS60))*2*self.pix_scale
        Ntot = np.array([Nx,Ny])
        mask = np.zeros([int(Ny),int(Nx)],dtype = bool)
        
        self.pix_coords = self.hex_centers*self.pix_scale + Ntot/2.
        # Lhex = np.array(self.hex_mask.shape)
        
        data = np.ones(len(self.valid_row),dtype=bool)
        
        for i in range(len(self.hex_centers)):
            # x_min = int(pix_coords[i,0] - Lhex[1]/2)
            # y_min = int(pix_coords[i,1] - Lhex[0]/2)
            # x_max = x_min + Lhex[1]
            # y_max = y_min + Lhex[0]
            # self.full_mask[y_min:y_max,x_min:x_max] *= self.hex_mask
            x_shift = self.pix_coords[i,0]
            y_shift = self.pix_coords[i,1]
            row = self.valid_row + y_shift
            col = self.valid_col + x_shift
            sparseMatrix = csr_matrix((data, (row,col)),  
                                      shape = mask.shape).toarray()
            # plt.figure()
            # plt.imshow(sparseMatrix,origin='lower')
            
            mask += sparseMatrix
            
        self.full_mask = ~mask
        
        
    def calculate_interaction_matrix(self,N_modes):
        """ Computes the interaction matrix: [n_pixels,n_hexes*n_modes] """

        self.int_mat = np.zeros([self.full_mask.size,len(self.hex_centers)*N_modes],dtype=np.uint8)
        
        L_x,L_y = np.shape(self.hex_mask)
        
        modes = np.zeros([L_x,L_y,N_modes])
        for j in range(N_modes):
            modes[:,:,j] = czern(j, [L_x,L_y])

        for i in range(len(self.hex_center_coords)):
            x_min = int(self.pix_coords[i,0] - L_y/2)
            y_min = int(self.pix_coords[i,1] - L_x/2)
            x_max = x_min + L_y
            y_max = y_min + L_x
            
            for j in range(N_modes):
                global_mode = np.zeros_like(self.full_mask)
                global_mode[y_min:y_max,x_min:x_max] = modes[:,:,j]
                self.int_mat[:,i*N_modes+j] = global_mode.flatten()


