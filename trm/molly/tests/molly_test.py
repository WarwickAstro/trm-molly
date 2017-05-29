import unittest
import copy
import tempfile

import numpy as np
from collections import OrderedDict as Odict
from trm import molly

class TestMolly(unittest.TestCase):
    """Provides simple tests of trm.molly classes and methods
    """

    def setUp(self):

        self.npix = 1500
        self.wave1 = 3800.
        self.wave2 = 5200.
        self.w = np.linspace(self.wave1,self.wave2,self.npix)
        self.f = np.ones_like(self.w)
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

        # make a spectrum
        self.spec = molly.Molly.fromdata(
            self.w,self.f,self.fe,self.cfrat,self.head
        )

    def test_fromdata(self):
        self.assertEqual(
            self.spec.narc, 2,
            'incorrect number of arc coefficients'
        )
        self.assertEqual(
            self.spec.npix, self.npix,
            'incorrect number of pixels'
        )
        self.assertAlmostEqual(
            self.spec.arc[0]+self.spec.arc[1]/self.npix, self.wave1,
            7, 'wavelength of first pixel wrong'
        )
        self.assertAlmostEqual(
            self.spec.arc[0]+self.spec.arc[1], self.wave2,
            7, 'wavelength of last pixel wrong'
        )

    def test_wave(self):
        w = self.spec.wave
        self.assertTrue((np.abs(w-self.w) < 1.e-4).all(),
                        'returned wavelengths are wrong')

    def test_write_read(self):
        # write to then read from a temporary file
        with tempfile.TemporaryFile() as fptr:
            self.spec.write(fptr)
            fptr.seek(0)
            spec = molly.Molly.fromfile(fptr)

        # check everything is recovered
        self.assertEqual(
            self.spec.npix, spec.npix, 'npix write/read problem'
        )
        self.assertEqual(
            self.spec.narc, spec.narc, 'narc write/read problem'
        )
        self.assertEqual(
            len(self.spec.head), len(spec.head), 'header write/read problem'
        )
        for key in self.spec.head.keys():
            self.assertTrue(key in spec.head,
                            'header key write/read problem')
            self.assertTrue(self.spec.head[key] == spec.head[key],
                            'header value write/read problem')
        self.assertTrue(
            (self.spec.wave == spec.wave).all(),
            'wavelengths write/read problem')
        self.assertTrue(
            (np.abs(self.spec.f - spec.f) < 1.e-6).all(),
            'fluxes write/read problem')

if __name__ == '__main__':
    unittest.main()
