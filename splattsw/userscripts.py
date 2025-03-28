from m4.mini_OTT import timehistory as th
import numpy as np


def analyze_opt_step(tn):
    fl = th.fileList(tn)
    cc = th.cubeFromList(fl)
    dd = cc-cc[0]
    zz = []
    for i in dd:
        zz.append(th.removeZernike(i,[1,2,3,4]))
    zz = np.ma.masked_array(zz)
    vv = zz.std( axis=(1,2))
    return vv

