import numpy as np
import matplotlib.pyplot as plt

from utilities import dm_system_setup
from utilities import fitting_error

config_tn = '20240920'
dm = dm_system_setup(config_tn)

seg0 = dm.segment[0]
local_fit_err = fitting_error(seg0.mask, seg0.IM, seg0.IFF, seg0.R)

global_fit_err = fitting_error(dm.global_mask, dm.glob_IM, dm.IFF, dm.R)
