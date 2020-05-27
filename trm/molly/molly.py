#!/usr/bin/env python

"""Package for reading and writing "molly" spectra.

"molly" is a program for the analysis of astronomical spectra which uses a
bespoke data format. This package reads molly spectra allowing access to
the internals of molly data. It can also write out spectra in molly format
if you need to import data into "molly".
"""

from collections import OrderedDict as Odict
import math, copy, struct
import numpy as np
from scipy import polyfit

__all__ = ['Molly', 'gmolly', 'rmolly', 'skip_molly', 'C']

# Speed of light in km/s
C = 299792.458

class Molly:
    """Container for a molly spectrum. Contains many somewhat obscure attributes
    closely aligned with their equivalents in the molly F77 program which
    should be consulted where the following is insufficient::

      fcode : int
           an integer format code. Molly spectra can differ as to whether
           there are errors and this code distinguishes different cases.

       head : OrderedDict
           ordered dictionary of header items. NB. molly's internal headers
           support a restricted variety of data types (4-byte integers, 4- and
           8-byte floats and strings (32 bytes at most). The names (key
           values) are at most 16 bytes. The limitations are enforced if you
           try to write to molly. Strings will be truncated and invalid data
           types ignored, so you should check that the headers inside molly
           are what you expect.

       npix : int
           number of pixels

       narc : int
           number of arc coefficients. Can be 0 meaning there are no arc
           coefficients. If positive, the coefficients are the direct poly
           coefficients representing the wavelength scale; if negative, they
           give the natural log of the wavelength scale.

       arc : numpy.array
          the arc polynomial coefficients (if narc ne 0). Wavelengths are
          always in Angstroms. If narc > 0 then the wavelengths are the sum
          over index i of arc[i]*(x**i) where x = n/npix with n the array
          element number starting with 1 (reflecting the FORTRAN legacy). If
          narc < 0, the exponential of this is taken. You should not normally
          need to know this if you use 'wave' to get the wavelengths, but you
          do if you want to manipulate the wavelength scale.

       label : bytes
           label to use to describe fluxes

       units : bytes
           the units of the fluxes. Should be 'MILLIJANSKYS', 'COUNTS' or
           'FLAMBDA'.

       f : numpy.ndarray
           fluxes (fcode == 3, 4 or 5) or counts(fcode = 1 or 2)

       fe : numpy.ndarray
           their uncertainties. < 0 indicates a point is masked. Can be
           None if there were no uncertainties.

       cfrat : numpy.ndarray
           the ratio of fluxes divided by counts. Can be None.

       wave  : numpy.ndarray
           heliocentric wavelengths.

    The best way to get a feel for these is to read in some molly spectra
    and look at these attributes one by one.
    """

    def __init__(self, fcode, head, npix, narc, arc,
                 label, units, f, fe, cfrat):
        """Initialises a Molly."""
        self.fcode = fcode
        self.head = head
        self.npix = npix
        self.narc = narc
        self.arc = arc
        self.label = label
        self.units = units
        self.f = f
        self.fe = fe
        self.cfrat = cfrat

    @classmethod
    def fromfile(cls, fptr):
        """Reads in a molly spectrum from a file-like object fptr, assumed
        to be set at the start of a spectrum. On successful execution, fptr
        will be at the start of the next spectrum.
        """

        fcode, head, npix, narc, arc, border = read_molly_head(fptr)
        label, units, f, fe, cfrat = read_molly_data(fptr, fcode, npix, border)

        return cls(fcode, head, npix, narc, arc,
                   label, units, f, fe, cfrat)

    @classmethod
    def fromdata(cls, w, f, fe=None, cfrat=None, head=None,
                 maxdev=0.001, nmax=6):
        """Creates a Molly from wavelength and flux arrays.

        Arguments::

          w : (numpy.array)
            wavelengths in Angstroms. Can be None. If the wavelengths are
            heliocentric then there should be no 'Vearth' header item or
            it should be set = 0.

          f : (numpy.array)
            array of fluxes.

          fe : (numpy.array)
            array of flux uncertainties

          cfrat : (numpy.array)
            ratios fluxes / counts

          head : (collections.OrderedDict)
            header names & values

          maxdev : (float)
            maximum deviation in terms of pixels when fitting wavelength
            scale.  Attempts will be made to fit polynomials to the wavelength
            scale.  First a linear wavelength scale, then a uniform
            logarithmic scale then a series of polynomials. If no fit ever has
            a maximum deviate below maxdev, an Exception will be raised.

          nmax : (int)
            maximum number of poly coefficients to try. e.g. If you wanted to
            stop at third-degree (cubic) polynomials, you would set nmax=4.

        """

        if fe is not None:
            fcode = 3
        else:
            fcode = 4

        npix = len(f)

        if w is None:
            narc = 0
            arc = np.array([])
        else:
            xpix = np.linspace(1,npix,npix)/float(npix)
            mdisp = (w.max()-w.min())/(npix - 1)

            # First attempt to fit the wavelengths with a linear fit
            arc = polyfit(xpix, w, 1)
            poly = np.poly1d(arc)
            dmax = abs(w-poly(xpix)).max()/mdisp
            dmin = dmax
            if dmax < maxdev:
                narc = 2
                arc = arc[::-1]
            else:
                # try a logarithmic scale
                lw = np.log(w)
                mldisp = (lw.max()-lw.min())/(npix - 1)
                arc = polyfit(xpix, lw, 1)
                poly = np.poly1d(arc)
                dmax = abs(lw-poly(xpix)).max()/mldisp
                if dmax < maxdev:
                    narc = -2
                    arc = arc[::-1]
                else:
                    dmin = min(dmin, dmax)
                    # linear fits have failed; try higher order polynomials
                    for ndeg in range(2,nmax):
                        arc = polyfit(xpix, w, ndeg)
                        poly = np.poly1d(arc)
                        dmax = abs(w-poly(xpix)).max()/mdisp
                        dmin = min(dmin, dmax)
                        if dmax < maxdev:
                            narc = ndeg + 1
                            arc = arc[::-1]
                            break
                    else:
                        raise Exception('No acceptable poly fit found; best case deviated by {0:f} pixels'.format(dmin))


        return cls(fcode, head, npix, narc, arc,
                   'fnu', 'MILLIJANSKYS', f, fe, cfrat)

    @property
    def wave(self):
        """Returns a (heliocentric) wavelength array in Angstroms.

        Raises an Exception is no scale is set. If a header item 'Vearth' is
        set, this is assumed to be the velocity (in km/s) the target acquires
        due to the motion of Earth relative to the Sun, and it will be removed
        from the wavelengths by multiplying the wavelengths by 1-v/c.
        """
        if self.narc != 0:
            x = np.polyval(self.arc[::-1], np.linspace(1,self.npix,self.npix)/self.npix)
            if self.narc < 0:
                x = np.exp(x)
            # correct to a heliocentric scale
            if 'Vearth' in self.head:
                x *= 1.-self.head['Vearth']/C
            return x
        else:
            raise Exception('No wavelength scale set')

    def write(self, fptr):
        """
        Writes a Molly to a file in molly format. The file must have been
        opened in binary mode and positioned at its end. At the end of running
        this routine, another Molly could be written to the same file to build
        up a large number of spectra.

        Arguments::

          fptr : (file-like object)
             Open as e.g. fptr = open('spec.mol','wb')

        """

        # a few checks because if we write invalid data it will be very hard
        # to diagnose.
        if self.narc != 0 and len(self.arc) != abs(self.narc):
            raise Exception('narc = {:d} does not match number of coefficients in arc (={:d})'.format(self.narc,len(self.arc)))

        if len(self.f) != self.npix:
            raise Exception('npix = {:d} does not match number of fluxes (={:d})'.format(self.npix,len(self.f)))

        if self.fe is not None and len(self.fe) != self.npix:
            raise Exception('npix = {:d} does not match number of flux errors (={:d})'.format(self.npix,len(self.fe)))

        if self.cfrat is not None and len(self.cfrat) != self.npix:
            raise Exception('npix = {:d} does not match number of flux/count ratios (={:d})'.format(self.npix,len(self.cfrat)))

        # first record. F77 writes the number of bytes in each
        # record at the start and end of the record.
        nbytes = 44

        # split headers into 4 supported types

        if self.head is not None:
            # character strings
            chead = Odict()
            for k,v in self.head.items():
                if isinstance(v, str):
                    chead[k] = v

            # 8-byte floats
            dhead  = Odict()
            for k,v in self.head.items():
                if isinstance(v, float):
                    dhead[k] = v

            # 4-byte ints
            ihead  = Odict()
            for k,v in self.head.items():
                if isinstance(v, int):
                    ihead[k] = v

            # 4-byte floats
            fhead  = Odict()
            for k,v in self.head.items():
                if isinstance(v, np.float32):
                    fhead[k] = v

        else:
            chead, dhead, ihead, fhead = Odict(), Odict(), Odict(), Odict()

        # write the record
        fptr.write(
            struct.pack(
                '2i16s7i',nbytes,self.fcode,tobytes(self.units,16),
                self.npix,self.narc,
                len(chead),len(dhead),len(ihead),len(fhead),
                nbytes)
        )

        # second record (names of header items)
        nbytes = 16*(len(chead)+len(dhead)+len(ihead)+len(fhead))
        fptr.write(struct.pack('i',nbytes))
        for key in chead.keys():
            fptr.write(struct.pack('16s',tobytes(key,16)))
        for key in dhead.keys():
            fptr.write(struct.pack('16s',tobytes(key,16)))
        for key in ihead.keys():
            fptr.write(struct.pack('16s',tobytes(key,16)))
        for key in fhead.keys():
            fptr.write(struct.pack('16s',tobytes(key,16)))
        fptr.write(struct.pack('i',nbytes))

        # third record (values of header items)
        nbytes = 32*len(chead) + 8*len(dhead) + 4*len(ihead) + 4*len(fhead)
        fptr.write(struct.pack('i',nbytes))
        for value in chead.values():
            fptr.write(struct.pack('32s',tobytes(value,32)))
        for value in dhead.values():
            fptr.write(struct.pack('d',value))
        for value in ihead.values():
            fptr.write(struct.pack('i',value))
        for value in fhead.values():
            fptr.write(struct.pack('f',value))
        fptr.write(struct.pack('i',nbytes))

        # fourth record (the arc coefficients)
        nbytes = 8*abs(self.narc)
        fptr.write(struct.pack('i',nbytes))
        if self.narc != 0:
            self.arc.tofile(fptr)
        fptr.write(struct.pack('i',nbytes))

        # fifth record (data)
        if self.fcode == 1 or self.fcode == 4:
            nbytes = 4*self.npix
            fptr.write(struct.pack('i',nbytes))
            np.cast['float32'](self.f).tofile(fptr)
            fptr.write(struct.pack('i',nbytes))

        elif self.fcode == 2 or self.fcode == 5:
            nbytes = 8*self.npix
            fptr.write(struct.pack('i',nbytes))
            np.cast['float32'](self.f).tofile(fptr)
            np.cast['float32'](self.fe).tofile(fptr)
            fptr.write(struct.pack('i',nbytes))

        elif self.fcode == 3:
            nbytes = 12*self.npix
            counts = self.f * self.cfrat
            errors = self.fe * self.cfrat
            flux = self.f.copy()
            print(np.linspace(1,len(flux),len(flux))[counts == 0.])

            flux[counts == 0.] = self.cfrat[counts == 0.]
            ff = np.cast['float32'](counts)
            print(
                flux[counts == 0.], counts[ff == 0.]
            )
            fptr.write(struct.pack('i',nbytes))
            np.cast['float32'](counts).tofile(fptr)
            np.cast['float32'](errors).tofile(fptr)
            np.cast['float32'](flux).tofile(fptr)
            fptr.write(struct.pack('i',nbytes))

        else:
            raise Exception('fcode = ' + str(self.fcode) + ' not implemented')

    def __len__(self):
        """Number of pixels in the Molly"""
        return self.npix

    def __repr__(self):
        return 'Molly(fcode={:d}, head={!r}, npix={:d}, narc={:d}, arc={!r}, label={:s}, units={:s}, f={!r}, fe={!r}, cfrat={!r})'.format(
            self.fcode, self.head, self.npix, self.narc, self.arc,
            self.label, self.units, self.f, self.fe, self.cfrat
        )

def gmolly(fname):
    """Generator for reading a molly file.

    Use as

      for spec in gmolly(fname):
         ... do something with spec
    """
    with open(fname, 'rb') as fptr:
        try:
            while True:
                yield Molly.fromfile(fptr)
        except EOFError:
            pass

def rmolly(fname, nspec=0):
    """Reads a spectrum or spectra from a molly file.

    Arguments::

       fname : (string)
          name of molly file

       nspec : (int)
          spectrum number to read. 0 = the whole lot returned as
          a list.
    """
    if nspec == 0:
        mlist = []
        for mspec in gmolly(fname):
            mlist.append(mspec)
        return mlist
    else:
        with open(fname, 'rb') as fptr:
            for n in range(nspec-1):
                skip_molly(fptr)
            mspec = Molly.fromfile(fptr)
        return mspec

def skip_molly(fptr):
    """Skips a molly spectrum.

    Assumes that the file pointer (opened in binary mode) is positioned at the
    start of the spectrum.  Leaves the pointer at the start of the next
    spectrum (if there is one).  Returns False if there is no spectrum to
    skip; raise Exceptions if there are other problems.

    """

    # If 'nbytes' in the next line comes up blank, we have reached the end of
    # the file
    nbytes = fptr.read(4)
    if nbytes == '': return False

    # If it does not start with 44 in either big or little endian form,
    # something is wrong
    (nbyte,)  = struct.unpack('<i', nbytes)
    if nbyte != 44:
        (nbyte,)  = struct.unpack('>i', nbytes)
        if nbyte != 44:
            raise Exception(
                'skip_molly: not a molly spectrum: first 4 bytes = {:d} not 44'.format(nbyte))
        border = '>'
    else:
        border = '<'

    # Read first line with various format items
    try:
        fcode,units,npix,narc,nchar,ndoub,nint,nfloat = struct.unpack(
            border + 'i16s6i',fptr.read(44))
    except:
        raise Exception("skip_molly: failed to read first line of molly spectrum")

    # compute number of bytes to skip.
    nskip = 36 + 16*(nchar+ndoub+nint+nfloat) + \
            32*nchar + 8*ndoub + 4*nint + 4*nfloat + 8*abs(narc)

    if fcode == 1:
        nskip += 4*npix
    elif fcode == 2:
        nskip += 8*npix
    elif fcode == 3:
        nskip += 12*npix
    elif fcode == 5:
        nskip += 8*npix
    else:
        raise Exception('skip_molly: invalid FCODE in molly spectrum = {:d}'.format(fcode))

    # skip them
    fptr.seek(nskip,1)
    return True

def tostr(barray):
    return barray.decode('utf-8')

def tobytes(strng, length):
    """Returns a string as a bytes object of fixed length (truncated
    or padded with blanks if need be"""
    b = strng.encode('utf-8')[:length]
    return b + (length-len(b))*b' '

def read_molly_head(fptr):
    """Reads the headers and arc (if present) of a molly spectrum.

    Returns (fcode, head, narc, arc, border) where fcode is the molly format
    code needed for reading the data, head is an collections.OrderedDict of
    header values, narc is the number of arc coefficients, arc are the arc
    coefficients, and border defines the byte order (either '>' or
    '<'). If narc=0, then arc will be None.

    Raises an EOFError if called when at the end of the file.
    raises Exceptions is other problems are encountered.

    """

    # If 'nbytes' in the next line comes up blank, we have reached the end of
    # the file
    nbytes = fptr.read(4)
    if nbytes == b'':
        raise EOFError('Reached end of fptr.name')

    # If it does not start with 44 in either big or little endian form,
    # something is wrong
    (nbyte,) = struct.unpack('<i', nbytes)
    if nbyte != 44:
        (nbyte,)  = struct.unpack('>i', nbytes)
        if nbyte != 44:
            raise Exception('read_molly_header: not a molly spectrum: first 4 bytes = ' + str(nbyte) + ' not 44')
        border = '>'
    else:
        border = '<'

    # Read first line with various format items
    try:
        fcode,units,npix,narc,nchar,ndoub,nint,nfloat = \
            struct.unpack(border + 'i16s6i',fptr.read(44))
    except:
        raise Exception("Failed to read first line of molly spectrum")

    # skip bytes at end of first record and at start of second
    fptr.seek(8,1)

    # read names of string header items
    cnames = []
    for i in range(nchar):
        name = tostr(fptr.read(16).strip())
        cnames.append(name)

    # read names of double header items
    dnames = []
    for i in range(ndoub):
        name = tostr(fptr.read(16).strip())
        dnames.append(name)

    # read names of integer header items
    inames = []
    for i in range(nint):
        name = tostr(fptr.read(16).strip())
        inames.append(name)

    # read names of float header items
    fnames = []
    for i in range(nfloat):
        name = tostr(fptr.read(16).strip())
        fnames.append(name)

    # skip bytes at end of second record and at start of third
    fptr.seek(8,1)

    # create header
    head = Odict()

    for i in range(nchar):
        value = tostr(fptr.read(32).strip())
        head[cnames[i]] = value

    dvals = struct.unpack(border + str(ndoub) + 'd', fptr.read(8*ndoub))
    for i in range(ndoub):
        head[dnames[i]] = dvals[i]

    ivals = struct.unpack(border + str(nint) + 'i', fptr.read(4*nint))
    for i in range(nint):
        head[inames[i]] = ivals[i]

    fvals = struct.unpack(border + str(nfloat) + 'f', fptr.read(4*nfloat))
    for i in range(nfloat):
        head[fnames[i]] = np.float32(fvals[i])

    # skip bytes at end of third record and at start of fourth
    fptr.seek(8,1)

    if narc != 0:
        arc = np.fromfile(file=fptr, dtype=border + 'f8', count=abs(narc))
    else:
        arc = None

    # skip 4 bytes at end of headers
    fptr.seek(4,1)

    return (fcode, head, npix, narc, arc, border)

def read_molly_data(fptr, fcode, npix, border):
    """label, units, f, fe, cfrat = read_molly_data(fptr, fcode, npix, border)

    Arguments::

      fptr : (file object)
        a file object opened in binary mode and positioned at the start of the
        data.

      fcode : (int)
        molly format code. Should have come from the header read stage

      npix : (int)
        number of pixels, from the header read.

      border : (string)
        string defining the byte order. Either '>' or '<'

    Reads data of a molly spectrum, assuming the header and arc have been read
    Returns label, units, fluxes, flux errors, and counts/flux ratio, some of
    which could come back as None.
    """

    # skip 4 bytes at start
    fptr.seek(4,1)

    cfrat = None

    if fcode == 1:
        f = np.fromfile(file=fptr, dtype=border + 'f4', count=npix)
        fe = None
        flabel = 'Counts'
        funits = 'COUNTS'

    elif fcode == 2:
        f = np.fromfile(file=fptr, dtype=border + 'f4', count=npix)
        fe = np.fromfile(file=fptr, dtype=border + 'f4', count=npix)
        flabel = 'Counts'
        funits = 'COUNTS'

    elif fcode == 3:
        counts = np.fromfile(file=fptr, dtype=border + 'f4', count=npix)
        errors = np.fromfile(file=fptr, dtype=border + 'f4', count=npix)
        flux = np.fromfile(file=fptr, dtype=border + 'f4', count=npix)

        cfrat = np.empty(npix, dtype=border + 'f4')
        mod = counts == 0.
        cfrat[mod] = flux[mod]
        mod = counts != 0.
        cfrat[mod] = counts[mod] / flux[mod]

        fe = np.empty_like(errors)
        ok = cfrat > 0.
        fe[ok] = errors[ok] / cfrat[ok]
        fe[~ok] = -1.
        f = flux
        f[counts == 0.] = 0.

        flabel = 'fnu'
        funits = 'MILLIJANSKYS'

    elif fcode == 4:
        f = np.fromfile(file=fptr, dtype=border + 'f4', count=npix)
        fe = None
        flabel = 'fnu'
        funits = 'MILLIJANSKYS'

    elif fcode == 5:
        f = np.fromfile(file=fptr, dtype=border + 'f4', count=npix)
        fe = np.fromfile(file=fptr, dtype=border + 'f4', count=npix)
        flabel = 'fnu'
        funits = 'MILLIJANSKYS'

    else:
        raise Exception(
            '_read_molly_data: invalid FCODE in molly spectrum = ' + str(fcode))

    # skip 4 bytes at end
    fptr.seek(4,1)

    return (flabel, funits, f, fe, cfrat)
