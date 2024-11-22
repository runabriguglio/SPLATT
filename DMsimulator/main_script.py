import numpy as np
import matplotlib.pyplot as plt

import utilities as utils

config_tn = '20240920'
dm = utils.dm_system_setup(config_tn)

# Segment Scramble
utils.segment_scramble(dm)

seg0 = dm.segment[0]
local_fit_err = utils.fitting_error(seg0.mask, seg0.IM, seg0.IFF, seg0.R)

# global_fit_err = utils.fitting_error(dm.global_mask, dm.glob_IM, dm.IFF, dm.R)


