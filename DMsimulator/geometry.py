import numpy as np
import matplotlib.pyplot as plt

from read_config import readConfig
from zernike_polynomials import computeZernike as czern

# from scipy.sparse import csr_matrix
# from scipy.sparse import random_array 
import scipy.sparse as sp

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
        
        mask_ul = np.fromfunction(lambda i,j: j >= L/2. + i/SIN60*COS60, [max_y,max_x])
        
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
        mask_shape = np.shape(mask)
        
        self.pix_coords = self.hex_centers*self.pix_scale + Ntot/2.

        # data = np.ones(len(self.valid_row),dtype=bool)
        
        # self.hex_indices = np.zeros([len(self.hex_centers),len(self.valid_row)],dtype=int)
        
        # for i in range(len(self.hex_centers)):
        #     row = (self.valid_row + self.pix_coords[i,1]).astype(int)
        #     col = (self.valid_col + self.pix_coords[i,0]).astype(int)
        #     sparse_mat = sp.csr_matrix((data, (row,col)),  
        #                               mask_shape, dtype=bool).toarray()
        #     # plt.figure()
        #     # plt.imshow(sparseMatrix,origin='lower')
        #     # print([i,np.max(row)*np.shape(mask)[1]+np.max(col)>np.size(mask)])
            
        #     # Save valid indices
        #     self.hex_indices[i,:] = col + row*mask_shape[1]
            
        #     mask += sparse_mat
        n_hex = n_hexagons(self.n_rings)
        rep_valid_row = np.tile(self.valid_row,n_hex)
        rep_valid_col = np.tile(self.valid_col,n_hex)
        
        valid_len = len(self.valid_row)
        rep_pix_coords = np.repeat(self.pix_coords, valid_len, axis = 0)
        
        row = (rep_valid_row + rep_pix_coords[:,1]).astype(int)
        col = (rep_valid_col + rep_pix_coords[:,0]).astype(int)
        
        data = np.ones(n_hex*valid_len,dtype=bool)
        
        sparse_mat = sp.csr_matrix((data, (row,col)),  
                                  mask_shape, dtype=bool)
          
        mask += sparse_mat
            
        self.full_mask = ~mask
        
        # Save valid indices
        hex_idx = col + row*mask_shape[1]
        self.hex_indices = np.reshape(hex_idx,[n_hex,valid_len])
        
        
    def calculate_interaction_matrix(self,n_modes):
        """ Computes the interaction matrix: 
            [n_pixels,n_hexes*n_modes] """
        n_hex = n_hexagons(self.n_rings)
        hex_data_len = np.sum(1-self.hex_mask)
        row_modes = np.zeros([hex_data_len*n_modes]) #,dtype = np.uint8
        for j in range(n_modes):
            row_modes[hex_data_len*j:hex_data_len*(j+1)] = czern(j+1, self.hex_mask)
            
        print('Computing interaction matrix...')      
        row_indices = np.tile(self.hex_indices,n_modes)
        mode_indices = np.arange(int(n_modes*n_hex))
        
        row = row_indices.flatten()
        col = np.repeat(mode_indices,hex_data_len)
        
        data = np.tile(row_modes,n_hex)

        int_mat_shape = [np.size(self.full_mask),n_modes*n_hex]
        self.int_mat = sp.csr_matrix((data, (row,col)),  
                                  int_mat_shape)
        
        
    def segment_scramble(self):
        """ Defines an initial segment scramble
        returning a masked array image """
        n_hex = n_hexagons(self.n_rings)
        
        # Retrieve number of modes from the interaction matrix
        n_modes = int(np.shape(self.int_mat)[1]/n_hex)
        
        # Generate random mode coefficients
        mode_vec = np.random.randn(n_hex*n_modes)
        
        # Matrix product
        print('Performing matrix product ...')
        # flat_img = (self.int_mat.dot(mode_vec)).astype(np.float16)
        flat_img = self.int_mat*mode_vec
        
        # Reshape and mask image
        print('Plotting....')
        img = np.reshape(flat_img, np.shape(self.full_mask))
        img = np.ma.masked_array(img, self.full_mask)
        
        return img
            
        
        
        
        
        
        
        
