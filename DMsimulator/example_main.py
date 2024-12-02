import numpy as np
import matplotlib.pyplot as plt

import utilities as utils

config_tn = '20240920'
# config_tn = 'JWST'
dm = utils.dm_system_setup(config_tn)

# Segment Scramble
utils.segment_scramble(dm, apply_shape=True)

# # Fitting errors
# seg0 = dm.segment[0]
# local_fit_err = utils.fitting_error_plots(seg0.mask, seg0.ZM, seg0.IFF, seg0.R)

# global_fit_err = utils.fitting_error_plots(dm.mask, dm.glob_ZM, dm.IFF, dm.R, dm.valid_ids.flatten())


# Update coordinates
utils.update_act_coords_on_ring(dm, 0)
utils.update_act_coords_on_ring(dm, 2)

# Apply flat and show actuators
dm.apply_flat()
dm.get_position(1)
dm.surface()

# # CapSens matrix
# meas_gap = utils.capsens_measure(dm, 1)