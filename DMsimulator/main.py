import numpy as np
import matplotlib.pyplot as plt
from geometry import Hexagons
# from astropy.io import fits
# from scipy.sparse import csr_matrix 
# from zernike_polynomials import computeZernike as czern

config_tn = '20240920'
hexes = Hexagons(config_tn)

mask = hexes.full_mask
plt.figure()
plt.imshow(mask,origin='lower',cmap='gray')

n_modes = 4
hexes.calculate_interaction_matrix(n_modes)
img = hexes.segment_scramble()
plt.figure()
# plt.grid('on')
plt.imshow(img,origin='lower',cmap='hot')

hexes.draw_hex_outline()


