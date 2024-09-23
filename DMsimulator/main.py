# import numpy as np
import matplotlib.pyplot as plt
from geometry import Hexagons
# from astropy.io import fits
# from scipy.sparse import csr_matrix 
# from zernike_polynomials import computeZernike as czern

config_tn = '20240920'
hexes = Hexagons(config_tn)

# Show the segments mask
mask = hexes.full_mask
plt.figure()
plt.imshow(mask,origin='lower',cmap='gray')

# Interaction matrix and initial scramble
n_modes = 11
hexes.calculate_interaction_matrix(n_modes)
img = hexes.segment_scramble()
plt.figure()
# plt.grid('on')
plt.imshow(img,origin='lower',cmap='hot')

# # Draw hexes and inscribed circle
# hexes.draw_hex_outline()

# # Test Zernike modes
# astig1 = czern(5,hexes.hex_mask)
# astig2 = czern(6,hexes.hex_mask)
# coma1 = czern(7,hexes.hex_mask)
# coma2 = czern(8,hexes.hex_mask)
# tref1 = czern(9,hexes.hex_mask)
# tref2 = czern(10,hexes.hex_mask)
# spher = czern(11,hexes.hex_mask)




