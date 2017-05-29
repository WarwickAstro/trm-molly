import unittest
import copy

import numpy as np
from collections import OrderedDict as Odict
from trm import molly

class TestMolly(unittest.TestCase):
    """Provides simple tests of CCD methods and attributes.

    """

    def setUp(self):

        self.npix = 1500
        xpix = np.linspace(1,self.npix,self.npix) / npix
        xpix -= xpix.mean()
        self.w = 4500. + 1000*xpix + 40.*xpix**2
        self.f = np.zeros_like(self.w)
        self.f += 1
        self.f += 0.6*np.exp(-((self.w-4600.)/30.)**2/2.)
        self.cfrat = 1000.*np.ones_like(self.w)
        counts = self.f*self.cfrat
        errors = np.sqrt(counts + 100)
        self.fe = errors / self.cfrat

        # Create header with all supported types
        self.head = Odict()
        self.head['RJD'] = 2451234.56789
        self.head['Object'] = 'Fake spectrum'
        self.head['Expose'] = np.float32(120.)
        self.head['Record'] = 12345

    def test_fromdata(self):
        spec = molly.Molly.fromdata(
            self.w,self.f,self.fe,self.cfrat,self.head)
        self.assertEqual(spec.narc,3,
                         'incorrect number of arc coefficients')

if __name__ == '__main__':
    unittest.main()
