#!/usr/bin/env python
import doctest
import unittest
import numpy as np
#from splattsw import webDAQ
#from splattsw.devices.utility import *
from splattsw.devices import webDAQ as wd
from splattsw import timehistory as th

sfreq = 1652. #to be adjusted or read from webDAQ or wdd file
#tt = timehistory.TimeHist()

aa = wd.WebDAQ() 
aa.connect()
fname = '/home/labot/Desktop/Data/ftp/SPLATT_2021-10-20T16-53-47-153.wdd'
gvect = wd.openwdd(fname)
spe, f = th.spectrum(gvect, dt=1/sfreq)

jb = 'SPLATT-OptBench'
aa.start_job(jobname=jb)
aa.get_jobs_status()
aa.get_schedule_status()
aa.read_data_current_job()

class webDAQTest(unittest.TestCase):

    def setUp(self):
        self.shape= (128, 128)
        self.radius = 52
        self.cx = 60.5
        self.cy = 70.5
        self.inradius=15
        self.params2check = [self.cx-0.5, self.cy-0.5, self.radius]
        self.testMask1= CircularMask(self.shape, 
                                     maskRadius=self.radius, 
                                     maskCenter=(self.cx, self.cy))
        self.testMask2= AnnularMask(self.shape, 
                                    maskRadius=self.radius, 
                                    maskCenter=(self.cx, self.cy), 
                                    inRadius=self.inradius )

    def test_docstring(self):
        import arte.utils.shape_fitter as pf_module
        doctest.testmod(pf_module, raise_on_error=True)
        
    def test_ransac(self):
        ff = ShapeFitter(self.testMask1.asTransmissionValue())
        ff.fit(method='RANSAC')
        ff2 = ShapeFitter(self.testMask2.asTransmissionValue())
        ff2.fit(method='RANSAC')
        self.assertEqual(np.round(self.params2check,1).all(), np.round(ff.params(),1).all())
        self.assertEqual(np.round(self.params2check,1).all(), np.round(ff2.params(),1).all())
    
    def test_correlation(self):
        ff = ShapeFitter(self.testMask1.asTransmissionValue())
        ff.fit()
        ff2 = ShapeFitter(self.testMask2.asTransmissionValue())
        ff2.fit()
        self.assertEqual(np.round(self.params2check,1).all(), np.round(ff.params(),1).all())
        self.assertEqual(np.round(self.params2check,1).all(), np.round(ff2.params(),1).all())
    
    def test_exceptions(self):
    
        ff = ShapeFitter(self.testMask1.asTransmissionValue())
    
        with self.assertRaises(Exception):
            _ = ff.fit(shape2fit='ellipse')
    
        with self.assertRaises(ValueError):
            _ = ff.fit(method='new')


       

if __name__ == "__main__":
    unittest.main()
