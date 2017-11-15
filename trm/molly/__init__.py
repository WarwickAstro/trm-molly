#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
trm.molly is a package for accessing spectra generated and written using
the F77 package "molly".

Basic usage::

  >> from trm import molly
  >>
  >> for spec in molly.gmolly('spectra.mol'):
  >>     print(spec)

to read in molly spectra from a file one by one, or

  >> from trm import molly
  >>
  >> specs = molly.rmolly('spectra.mol')
  >> print(specs)

to read in all spectra into a list. The object 'spec' above
is a trm.molly.Molly. Look at that for more detail. Many of
its attributes can be ignored, but the wavelengths, fluxes
and their errors can be accessed through attributes 'wave', 'f' 
and 'fe'



"""

from .molly import *

__all__ = molly.__all__
