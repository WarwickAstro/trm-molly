import sys
import os
import numpy as np
import pandas as pd
import time
from trm import molly

def mkmollydb(args=None):
    """Starting from the present working directory, this descends the
    tree to search for molly files, reads them, and creates a database
    from them including the location, spectrum number, object name,
    position, wavelength range etc.

    """

    dfs = []
    for dpath, dnames, fnames in os.walk('.'):
        mnames = [mname for mname in fnames if mname.endswith('.mol')]
        for mname in mnames:
            # generate full path to the molly file from pwd
            fullpath = os.path.join(dpath, mname)

            # extract data from the file
            fullpaths, objects, npixs, narcs, ras, decs, wmins, wmaxs, wdisps = ([] for i in range(9))
            try:
                for mspec in molly.gmolly(fullpath):
                    fullpaths.append(fullpath)
                    objects.append(mspec.head.get('Object',None))
                    ras.append(mspec.head.get('RA',None))
                    decs.append(mspec.head.get('Dec',None))
                    npixs.append(mspec.npix)
                    narcs.append(mspec.narc)
                    if mspec.narc != 0:
                        wave = mspec.wave
                        wmins.append(wave.min())
                        wmaxs.append(wave.max())
                        wdisps.append((wave.max()-wave.min())/(mspec.npix-1))
                    else:
                        wmins.append(None)
                        wmaxs.append(None)
                        wdisps.append(None)

            except Exception as err:
                print('Exception=',err)
                print(fullpath)

            # generate spectrum numbers
            nspecs = np.linspace(1,len(objects),len(objects)).astype(int)

            # create data dictionary
            ddict = {
                'Path' : fullpaths,
                'Nspec' : nspecs,
                'Object' : objects,
                'RA' : ras,
                'Dec' : decs,
                'Npix' : np.array(npixs,dtype=int),
                'Narc' : np.array(narcs,dtype=int),
                'Wmin' : wmins,
                'Wmax' : wmaxs,
                'Wdisp' : wdisps,
            }

            # add into list of DataFrames
            dfs.append(pd.DataFrame(ddict))

    bdf = pd.concat(dfs)
    print(bdf)
    bdf.to_pickle('molly.db')
