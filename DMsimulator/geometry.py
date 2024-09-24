import numpy as np
import matplotlib.pyplot as plt # only used in draw_hex_outline
from scipy.sparse import csr_matrix
import os

from read_config import readConfig
from zernike_polynomials import computeZernike as czern
from rotate_coordinates import rotate_by_60deg as crot60

from astropy.io import fits
from read_and_write_fits import write_to_fits


# Useful variables
SIN60 = np.sin(np.pi/3.)
COS60 = np.cos(np.pi/3.)


def n_hexagons(n_rings):
    return int(1 + (6 + n_rings*6)*n_rings/2)

class Hexagons():

    def __init__(self, TN):

        # Read configuration files
        config_path = './config_par_' + TN + '.yaml'
        dm_par, opt_par = readConfig(config_path)

        self.gap = dm_par[0]
        self.hex_side_len = dm_par[1]
        self.n_rings = dm_par[2]

        self.pix_scale = opt_par[0]
        
        self.savepath = './' + TN + '/'
        
        try: # Create new folder
            os.mkdir(TN)
        except FileExistsError: # Folder already existing
            pass
        
        self.define_hex_centers()
        self.define_local_mask()
        self.define_hexagon_valid_indices()
        self.define_global_mask()
        
    

    def define_hex_centers(self):
        """ Defines and saves the coordinates of the 
        centers of all hexagonal segments """
        
        file_path = self.savepath + 'hex_centers_coords.fits'
        try:
            with fits.open(file_path) as hdu:
                self.hex_centers = np.array(hdu[0].data)
            return
        
        except FileNotFoundError:
            pass
        
        # Number of hexes
        self.hex_centers = np.zeros([n_hexagons(self.n_rings),2])
        
        # Height of hex + gap
        L = self.gap + 2.*self.hex_side_len*SIN60
        
        for ring_ctr in range(self.n_rings):
            
            hex_ctr = n_hexagons(ring_ctr)
            ring_ctr += 1
            R = L*ring_ctr
            
            aux = np.array([R*SIN60,R*COS60])
            self.hex_centers[hex_ctr,:] = aux
            
            for i in range(5):
                aux = crot60(aux)
                self.hex_centers[hex_ctr+(i+1)*ring_ctr,:] = aux
            
            if ring_ctr > 1:
                for j in range(ring_ctr-1):
                    shift = self.gap + 2.*self.hex_side_len*SIN60
                    aux[0] = self.hex_centers[hex_ctr,0] 
                    aux[1] = self.hex_centers[hex_ctr,1] - (j+1)*shift
                    self.hex_centers[hex_ctr+j+1,:] = aux
                    for i in range(5):
                        aux = crot60(aux)
                        self.hex_centers[hex_ctr+j+1+(i+1)*ring_ctr,:] = aux

        # Save to fits
        write_to_fits(self.hex_centers,file_path)
        

    def define_local_mask(self):
        """ Forms the hexagonal mask to be placed 
        at the given hexagonal coordinates """
        
        file_path = self.savepath + 'hexagon_mask.fits'
        try:
            with fits.open(file_path) as hdu:
                self.local_mask = np.array(hdu[0].data).astype(bool)
            return
        
        except FileNotFoundError:
            pass
        
        L = self.hex_side_len*self.pix_scale
        
        # Consider points in the upper-left of the hexagon
        max_y = int(L*SIN60)
        max_x = int(L/2.+L*COS60)
        
        mask_ul = np.fromfunction(lambda i,j: j >= L/2. + i/SIN60*COS60, [max_y,max_x])
        
        mask = np.zeros([2*max_y,2*max_x], dtype=bool)
        
        mask[0:max_y,max_x:] = mask_ul # upper left
        mask[0:max_y,0:max_x] = np.flip(mask_ul,1) # upper right
        mask[max_y:,:] = np.flip(mask[0:max_y,:],0) # lower
                    
        self.local_mask = mask
        
        # Save to fits
        write_to_fits((self.local_mask).astype(np.uint8),file_path)
        
        
        
    def define_hexagon_valid_indices(self):
        """ Finds the full aperture image (containing all segments)
        row and column indices for the hexagons data """
        
        file_path = self.savepath + 'segments_indices.fits'
        try:
            with fits.open(file_path) as hdu:
                self.hex_indices = np.array(hdu[0].data)
                self.global_row_idx = np.array(hdu[1].data)
                self.global_col_idx = np.array(hdu[2].data)
            return
        
        except FileNotFoundError:
            pass
        
        # Height of hex + gap
        L = self.gap + 2.*self.hex_side_len*SIN60
        
        # Full mask dimensions
        Ny = (L*self.n_rings + self.hex_side_len*SIN60)*2*self.pix_scale
        Nx = (L*self.n_rings*SIN60 + self.hex_side_len*(0.5+COS60))*2*self.pix_scale
        Ntot = np.array([Nx,Ny])

        # Hexagon centers pixel coordinates
        self.pix_coords = self.hex_centers*self.pix_scale + Ntot/2.
        
        My,Mx = np.shape(self.local_mask)
        x = np.arange(Mx,dtype=int)
        y = np.arange(My,dtype=int)
        X,Y = np.meshgrid(x,y)
        local_X = X[~self.local_mask]
        local_Y = Y[~self.local_mask]
        local_row_idx = local_Y - int(My/2)
        local_col_idx = local_X - int(Mx/2)

        n_hex = n_hexagons(self.n_rings)
        rep_local_row = np.tile(local_row_idx,n_hex)
        rep_local_col = np.tile(local_col_idx,n_hex)
        
        # valid_len = len(local_row_idx)
        valid_len = np.sum(1-self.local_mask)
        rep_pix_coords = np.repeat(self.pix_coords, valid_len, axis = 0)
        
        self.global_row_idx = (rep_local_row + rep_pix_coords[:,1]).astype(int)
        self.global_col_idx = (rep_local_col + rep_pix_coords[:,0]).astype(int)
        
        # Save valid hexagon indices
        hex_idx = self.global_col_idx + self.global_row_idx*int(Nx)
        self.hex_indices = np.reshape(hex_idx,[n_hex,valid_len])
        
        # Save to fits
        write_to_fits([self.hex_indices,self.global_row_idx,self.global_col_idx],file_path)
        
        # # Save in a matrix, with dimensions [n_hex,2,valid_len]
        # hex_mat_idx = np.reshape(self.hex_indices,[2,n_hex,valid_len])
        # hex_mat_idx = np.transpose(hex_mat_idx,[1,0,2])
        
        # return hex_mat_idx

        
    def define_global_mask(self):
        """ Forms the segment mask """
        
        file_path = self.savepath + 'segments_mask.fits'
        try:
            with fits.open(file_path) as hdu:
                self.global_mask = np.array(hdu[0].data).astype(bool)
            return
        
        except FileNotFoundError:
            pass
        

        # Height of hex + gap
        L = self.gap + 2.*self.hex_side_len*SIN60
        
        # Full mask dimensions
        Ny = (L*self.n_rings + self.hex_side_len*SIN60)*2*self.pix_scale
        Nx = (L*self.n_rings*SIN60 + self.hex_side_len*(0.5+COS60))*2*self.pix_scale
        mask_shape = np.array([int(Ny),int(Nx)])
        mask = np.zeros(mask_shape,dtype = bool)
        
        #valid_len = len(self.global_row_idx)
        valid_len = np.sum(1-self.global_mask)
        
        data = np.ones(valid_len, dtype=bool)
        
        sparse_mat = csr_matrix((data, (self.global_row_idx,self.global_col_idx)),  
                                  mask_shape, dtype=bool)
        
        mask += sparse_mat
            
        self.global_mask = ~mask
        
        # Save to fits
        write_to_fits((self.global_mask).astype(np.uint8),file_path)
        
        
        
    def calculate_interaction_matrix(self,n_modes):
        """ Computes the interaction matrix: 
            [n_pixels,n_hexes*n_modes] """
            
        n_hex = n_hexagons(self.n_rings)
        int_mat_shape = [np.size(self.global_mask),n_modes*n_hex]
            
        file_path = self.savepath + str(n_modes) + 'modes_interaction_matrix.fits'
        try:
            with fits.open(file_path) as hdu:
                # self.int_mat = csr_matrix(hdu[0].data)
                mat_data = hdu[0].data
                indices = hdu[1].data
                indptr = hdu[2].data
                self.int_mat = csr_matrix((mat_data,indices,indptr), int_mat_shape)
            return
        
        except FileNotFoundError:
            pass
        
        hex_data_len = np.sum(1-self.local_mask)
        row_modes = np.zeros([hex_data_len*n_modes]) #,dtype = np.uint8
        for j in range(n_modes):
            row_modes[hex_data_len*j:hex_data_len*(j+1)] = czern(j+1, self.local_mask)
            
        print('Computing interaction matrix...')      
        row_indices = np.tile(self.hex_indices,n_modes)
        mode_indices = np.arange(int(n_modes*n_hex))
        
        row = row_indices.flatten()
        col = np.repeat(mode_indices,hex_data_len)
        
        data = np.tile(row_modes,n_hex)

        self.int_mat = csr_matrix((data, (row,col)),  
                                  int_mat_shape)
        
        # Save to fits
        print('Saving interaction matrix...') 
        data_list = []
        data_list.append(self.int_mat.data)
        data_list.append(self.int_mat.indices)
        data_list.append(self.int_mat.indptr)
        write_to_fits(data_list, file_path)
        
        
        
        
    def segment_scramble(self):
        """ Defines an initial segment scramble
        returning a masked array image """
        
        file_path = self.savepath + 'segment_scramble.fits'
        try:
            with fits.open(file_path) as hdu:
                img = np.array(hdu[0].data)
                img_mask = np.array(hdu[1].data).astype(bool)
            masked_img = np.ma.masked_array(img,mask=img_mask)
            return masked_img
        except FileNotFoundError:
            pass
        
        n_hex = n_hexagons(self.n_rings)
        
        # Retrieve number of modes from the interaction matrix
        n_modes = int(np.shape(self.int_mat)[1]/n_hex)
        
        # Generate random mode coefficients
        mode_vec = np.random.randn(n_hex*n_modes)
        
        # Probability inversely proportional to spatial frequency
        m = int(np.round((np.sqrt(8*n_modes)-1.)/2.))
        freq_vec = np.repeat(np.arange(m)+1,np.arange(m)+1)
        prob_vec = 1./freq_vec[0:n_modes]
        prob_vec_rep = np.tile(prob_vec,n_hex)
        
        # Modulate on the probability
        mode_vec = mode_vec * prob_vec_rep
        
        # Matrix product
        print('Performing matrix product...')
        flat_img = self.int_mat*mode_vec
        
        # Reshape and mask image
        print('Plotting....')
        img = np.reshape(flat_img, np.shape(self.global_mask))
        masked_img = np.ma.masked_array(img, self.global_mask)
        
        # Save to fits
        write_to_fits(masked_img, file_path)
        
        return masked_img
    
    
    def draw_hex_outline(self):
        """ Plots the hexagons' outline and the inscribed circle """
        
        hex_sides = np.zeros([8,2])
        hex_sides[0,:] = np.array([-2*COS60, 0.])
        hex_sides[1,:] = np.array([-0.5, SIN60])
        hex_sides[2,:] = np.array([-hex_sides[1,0],hex_sides[1,1]])
        hex_sides[3,:] = -hex_sides[0,:]
        hex_sides[4,:] = -hex_sides[1,:]
        hex_sides[5,:] = -hex_sides[2,:]
        hex_sides[6,:] = hex_sides[0,:]
        hex_sides[-1,:] = np.array([None,None])
        
        plt.figure()
        plt.grid('on')
        
        rep_c_coords = np.tile(self.hex_centers,len(hex_sides))
        rep_c_coords = rep_c_coords.flatten()
        hex_sides = hex_sides.flatten()
        rep_hex_sides = np.tile(hex_sides,len(self.hex_centers))
        coords = rep_c_coords + rep_hex_sides 
        coords = np.reshape(coords,[int(len(coords)/2),2])
        plt.plot(coords[:,0],coords[:,1],color='goldenrod')
                 
        
        # Plot inscribed cirle
        L = self.gap + 2.*self.hex_side_len*SIN60
        R = np.sqrt((L*self.n_rings)**2 + (self.hex_side_len*(0.5+COS60))**2) - self.hex_side_len*COS60
        x_vec = np.linspace(-R,R,100)
        y_vec = np.sqrt(R**2-x_vec**2)
        
        plt.plot(x_vec,y_vec,'--',color='green')
        plt.plot(x_vec,-y_vec,'--',color='green')
            
        plt.axis('equal')
        
        return coords

            
        
        
        
        
        
        
        
