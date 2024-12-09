import numpy as np
from M4.m4.mOTT_analysis import timehistory as th

basepath = '/home/labot/testfold/'
tn='tn0'
z2fit = np.array([1,2,3])
freq4d = 90.

flist=th.fileList(tn, fold=basepath)
zvec = th.zernikePlot(fl, modes=z2fit)
tvec = arange(len(flist))/freq4d





#def timevec(filelist):
#    np = len(filelist)
#    t = arange(np)/freq4d
#    return t

