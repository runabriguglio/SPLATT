# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 09:32:38 2024

@author: menes
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy import ndimage # For image rotation
# import cv2

# Define hexagon of unit side coordinates, exploiting symmetries
hex_sides = np.zeros([7,2])

COS60 = np.cos(np.pi/3.)
SIN60 = np.sin(np.pi/3.)
TAN60 = SIN60/COS60

hex_sides[0,:] = np.array([-2*COS60, 0.])
hex_sides[1,:] = np.array([-0.5, SIN60])
hex_sides[2,:] = np.array([-hex_sides[1,0],hex_sides[1,1]])
hex_sides[3,:] = -hex_sides[0,:]
hex_sides[4,:] = -hex_sides[1,:]
hex_sides[5,:] = -hex_sides[2,:]
hex_sides[-1,:] = hex_sides[0,:] # only used for plotting


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


class Mask:
    
    def __init__(self, number_of_rings, hex_size, gap, pixel_scale, cw_rot_angle = 0., center_point = np.array([0.,0.])):
        
        self.n_rings = number_of_rings
        self.hex_side_len = hex_size
        self.g = gap
        self.pix_scale = pixel_scale
        
        self.center = center_point
        self.angle = cw_rot_angle
        
        self._define_pixel_mask()
        
        
    def define_hex_centers(self):
        # Define the hexagons' center points
        ring_ctr = 0
        origin = np.array([0.,0.])
        
        hex_center_list = np.zeros([1,2])
        
        L = self.g + 2.*self.hex_side_len*SIN60
        
        for ring_ctr in range(self.n_rings):
            # tot_len += L
            ring_ctr += 1
            n_hexes = 6*ring_ctr
            hex_centers = np.zeros([n_hexes,2])
            
            R = L*ring_ctr
            
            aux = origin.copy()
            
            aux[0] += R*SIN60
            aux[1] += R*COS60
            hex_centers[0,:] = aux
            
            for i in range(5):
                aux = rotate_by_60deg(aux)
                hex_centers[(i+1)*ring_ctr,:] = aux
            
            if ring_ctr > 1:
                for j in range(ring_ctr-1):
                    shift = self.g + 2.*self.hex_side_len*SIN60
                    aux[0] = hex_centers[0,0] 
                    aux[1] = hex_centers[0,1] - (j+1)*shift
                    hex_centers[j+1,:] = aux.transpose()
                    for i in range(5):
                        aux = rotate_by_60deg(aux)
                        hex_centers[j+1+(i+1)*ring_ctr,:] = aux
                
            # hex_center_list.append(hex_centers)
            hex_center_list = np.concatenate([hex_center_list,hex_centers])
            
        self.hex_center_coords = np.array(hex_center_list)
            
        # Define full mask
        N_y = int((L*self.n_rings + self.hex_side_len*SIN60)*2*self.pix_scale)
        N_x = int((L*self.n_rings*SIN60 + self.hex_side_len*(0.5+COS60))*2*self.pix_scale)
        
        self.N = np.array([N_x,N_y])
        self.full_mask = np.zeros([self.N[1],self.N[0]])#,np.uint8)
        
        # Inscribed circle
        Radius = np.sqrt((L*self.n_rings)**2 + (self.hex_side_len*(0.5+COS60))**2) - self.hex_side_len*COS60
        
        # Circular Mask
        img_size = self.full_mask.shape
        xc = int(img_size[1]/2.)
        yc = int(img_size[0]/2.)
        r = Radius * self.pix_scale
        self.circle_mask = np.fromfunction(lambda i,j: (i-yc)**2 + (j-xc)**2 <= r**2, [self.N[1],self.N[0]])
    
    
    def draw_mask(self):
        # Draw the hexagons' perimeters
        plt.figure()
        plt.grid('on')
        
        for i in range(len(self.hex_center_coords)):
            c_coords = self.hex_center_coords[i] + self.center
            coords = c_coords + self.hex_side_len*hex_sides
            coords = cw_rotate(coords, self.angle)
            plt.plot(coords[:,0],coords[:,1],color='goldenrod')
        
        # Plot inscribed cirle
        L = self.g + 2.*self.hex_side_len*SIN60
        R = np.sqrt((L*self.n_rings)**2 + (self.hex_side_len*(0.5+COS60))**2) - self.hex_side_len*COS60
        x_vec = np.linspace(-R,R,100)
        y_vec = np.sqrt(R**2-x_vec**2)
        
        plt.plot(x_vec,y_vec,'--',color='green')
        plt.plot(x_vec,-y_vec,'--',color='green')
            
        plt.axis('equal')
        plt.show()
        
        plt.figure()
        plt.imshow(self.mask,origin='lower')
        
        
    def mask_pixels(self):
        
        L_x = np.shape(self.pixel_mask)[0]
        L_y = np.shape(self.pixel_mask)[1]
       
        c_coords = (self.hex_center_coords + self.center)*self.pix_scale + self.N/2.
        for i in range(len(self.hex_center_coords)):
            x_min = int(c_coords[i,0] - L_y/2)
            y_min = int(c_coords[i,1] - L_x/2)
            x_max = x_min + L_y
            y_max = y_min + L_x
            self.full_mask[y_min:y_max,x_min:x_max] += self.pixel_mask
            
        self.mask = np.multiply(self.full_mask,self.circle_mask)
        
        if np.abs(self.angle) > 0.:
            mask_size = np.array(self.mask.shape)
            # print(mask_size)
            rot_mask = ndimage.rotate(self.mask, self.angle/np.pi*180.)
            new_size = np.array(rot_mask.shape)
            shape_diff = new_size - mask_size
            shift = np.array([int(shape_diff[0]/2.),int(shape_diff[1]/2.)])
            self.mask = rot_mask[shift[0]:-shift[0],shift[1]:-shift[1]]
        
        
        
    def _define_pixel_mask(self):
        L = self.hex_side_len*self.pix_scale
        
        # Include points in the upper-left of the hexagon
        max_y = int(L*SIN60)
        max_x = int(L/2.+L*COS60)
        
        mask_ul = np.fromfunction(lambda i,j: j <= np.floor(L/2. + i/TAN60), [max_y,max_x])
        
        mask = np.zeros([2*max_y,2*max_x], np.uint8)
        
        mask[0:max_y,max_x:] = mask_ul # upper left
        mask[0:max_y,0:max_x] = np.flip(mask_ul,1) # upper right
        mask[max_y:,:] = np.flip(mask[0:max_y,:],0) # lower
                    
        self.pixel_mask = mask
        
        
    def project_Zernike(self, Z_number):
        mask_size = np.array(self.pixel_mask.shape)
        Y = mask_size[1]
        X = mask_size[0]
        
        # Normalization TBD: zero mean and std = 1
        match Z_number:
            case 0: # Piston
                mode = np.ones([X,Y])
            case 1: # Tip
                mode = np.fromfunction(lambda i,j: 2.*j/Y - 1., [X,Y])
            case 2: # Tilt
                mode = np.fromfunction(lambda i,j: 2.*i/X - 1., [X,Y])
            case 3: # Focus
                mode = np.fromfunction(lambda i,j: ((j-Y/2.)/Y)**2+((i-X/2.)/Y)**2, [X,Y])
            case _:
                print("Zernike mode too high: not yet implemented")
                return
            
        masked_mode = np.multiply(self.pixel_mask, mode)
        
        # plt.figure()
        # plt.imshow(masked_mode,origin='lower')
            
        return masked_mode
    
    
    def interaction_matrix(self, N_modes):   
        
        self.int_mat = np.zeros([self.full_mask.size,len(self.hex_center_coords)*N_modes])
        
        L_x = np.shape(self.pixel_mask)[0]
        L_y = np.shape(self.pixel_mask)[1]
        
        modes = np.zeros([L_x,L_y,N_modes])
        for j in range(N_modes):
            modes[:,:,j] = self.project_Zernike(j)
       
        c_coords = (self.hex_center_coords + self.center)*self.pix_scale + self.N/2.
        for i in range(len(self.hex_center_coords)):
            x_min = int(c_coords[i,0] - L_y/2)
            y_min = int(c_coords[i,1] - L_x/2)
            x_max = x_min + L_y
            y_max = y_min + L_x
            
            for j in range(N_modes):
                global_mode = np.zeros_like(self.full_mask)
                global_mode[y_min:y_max,x_min:x_max] = modes[:,:,j]
                self.int_mat[:,i*N_modes+j] = global_mode.flatten()
                
        
            
                
            


pixel_scale = 100.    # e.g. pixel/m
hex_side = 2.  
gap = 0.4
n_rings = 3
cw_rot_angle = 0

mask1 = Mask(n_rings, hex_side, gap, pixel_scale, cw_rot_angle)
mask1.define_hex_centers()
mask1.mask_pixels()
# mask1.draw_mask()
mask1.interaction_matrix(3)

# masked_ar = np.ma.masked_array(data,mask)

n_hex = len(mask1.hex_center_coords)
mode_vec = np.zeros([n_hex*3])

for i in range(n_hex):
    for j in range(3):
        mode_vec[i*3+j] = np.random.randn()
        
flat_img = mask1.int_mat @ mode_vec

img = np.reshape(flat_img, mask1.full_mask.shape)

plt.imshow(img)
        
        
