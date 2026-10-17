"""
Generate .FITS data files containing magnitude standard deviation, mean magnitude, and mean magnitude error from master
catalog tables.

Called from plot.py in /src/ directory.

Usage:
~/octarine/src$ python plot.py --std
"""

import os
import numpy as np
import re
from astropy.table import Table, Column

from daomop import storage

MAG_STD_DATA_DIRECTORY = 'plotting/mag_std_data/'


def std(hpx):
    """
    Generate a table containing the magnitude standard deviation, mean magnitude value, and mean magnitude error value.
    Puts the resulting table to VOSpace in the from of a FITS file.

    :param hpx: catalog HEALPIX
    """
    mag_filename = MAG_STD_DATA_DIRECTORY + str(hpx) + '_magnitudes.fits'
    std_filename = MAG_STD_DATA_DIRECTORY + str(hpx) + '_stds.fits'

    if not os.path.exists(std_filename) or not os.path.exists(mag_filename):
        cat = storage.HPXCatalog(hpx, catalog_dir='catalogs/master', dest_directory='master')
        table = cat.table
        art = storage.Artifact(storage.Observation(int(hpx)), version="", subdir='catalogs', ext='stds.fits')

        condition = (table['MATCHES'] > 2)
        hpxids = np.unique(table['HPXID'][condition])
        print len(hpxids)
        mags = np.array([table[table['HPXID'] == hpxid]['MAG_AUTO'] for hpxid in hpxids])
        mag_err = np.array([table[table['HPXID'] == hpxid]['MAGERR_AUTO'] for hpxid in hpxids])

        mean_magerr = [(np.mean(x)) for x in mag_err]
        stds = [(np.std(x)) for x in mags]
        mean_mag = [(np.mean(x)) for x in mags]

        table = Table()
        table.add_column(Column(stds, name='mag_std'))
        table.add_column(Column(mean_mag, name='mag_mean'))
        table.add_column(Column(mean_magerr, name='magerr_mean'))
        table.write(art.filename, format='fits', overwrite=True)
        art.put()


def main():
    if not os.path.exists(MAG_STD_DATA_DIRECTORY):
        os.mkdir(MAG_STD_DATA_DIRECTORY)

    directory = storage.listdir(os.path.join(os.path.dirname(storage.DBIMAGES),
                                             storage.CATALOG, 'master'), force=True)
    reg = []
    for item in directory:
        # HPX_02434_RA_185.6_DEC_+37.2
        x = re.match('HPX_(?P<pix>\d{5})_RA_(?P<ra>\d{3}\.\d)_DEC_\+(?P<dec>\d{2}\.\d)', item)

        hpx = int(x.group('pix'))
        try:
            if hpx not in reg and hpx > 0:  # alter 'hpx > 0' to use specific files
                reg.append(hpx)
                std(hpx)
        except TypeError as err:
            print err
