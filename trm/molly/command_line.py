import sys
import os
import math
from collections import OrderedDict as Odict

import numpy as np
from trm import molly

def genspec(args=None):
    """
    Generates fake spectra and save them to a molly file
    """

    # Create wavelength scale and flux/count ratio
    npix = 1500
    xpix = np.linspace(1,npix,npix) / npix
    xpix -= xpix.mean()
    w = 4500. + 1000*xpix + 40.*xpix**2
    cfrat = 1000.*np.ones_like(w)

    # Create header with all supported types
    head = Odict()
    head['Object'] = 'Fake spectrum'
    head['Expose'] = np.float32(120.)
    head['Record'] = 12345

    # Generate the trail
    NSPEC = 100
    P = 0.2
    K = 250.
    T0 = 2451234.56789
    STEP = 0.002
    with open('genspec.mol','wb') as fout:
        for n in range(NSPEC):
            t = head['RJD'] = T0 + STEP*n
            v = K*math.sin(2.*math.pi*(t-T0)/P)
            f = np.ones_like(w)
            f += 0.6*np.exp(-((w-4600.*(1+v/molly.C))/30.)**2/2.)
            counts = f*cfrat
            errors = np.sqrt(counts + 100)
            fe = errors / self.cfrat
            f = np.random.normal(f,fe)
            mspec = molly.Molly.fromdata(w,f,fe,cfrat,head)
            mspec.write(fout)
