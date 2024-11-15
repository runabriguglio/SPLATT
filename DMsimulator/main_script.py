import numpy as np
import matplotlib.pyplot as plt

from utilities import dm_system_setup
from utilities import fitting_error

config_tn = '20240920'
dm = dm_system_setup(config_tn)

seg0 = dm.segment[0]
local_fit_err = fitting_error(seg0.mask, seg0.IM, seg0.IFF, seg0.R)

rec_IFF = (dm.IFF[dm.valid_ids[0],:37]).todense()
rec_R = (dm.R[:37,dm.valid_ids[0]]).todense()

IFF_diff = seg0.IFF - rec_IFF
print(np.max(IFF_diff))
plt.figure()
plt.imshow(rec_R*rec_IFF)
plt.colorbar()

global_fit_err = fitting_error(dm.global_mask, dm.glob_IM, dm.IFF, dm.R)
