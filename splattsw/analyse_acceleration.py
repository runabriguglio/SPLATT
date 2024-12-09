from SPLATT.splattsw import splatt_analysis as sp

import matplotlib.pyplot as plt
import numpy as np

sp.wdsync()

wdfile = sp.last_wdfile()
data = sp.openfile(wdfile)
sp.plot_data(data)