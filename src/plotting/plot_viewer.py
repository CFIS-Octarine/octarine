"""
General image plotting/visualization script. Called from plot.py using optional arguments to specify which plots to
create. Depends on data files being stored locally in correct directories.
"""

import os
import sys
import re
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import wcs

from daomop import storage

HELIO_DATA_DIR = 'plotting/heliocentric_data/'
FITS_DATA_DIR = 'plotting/fits_data/'
MAG_STD_DIR = 'plotting/mag_std_data/'


def load_images(hpx, hpx_files):
    """
    Generate a figure consisting of 3 plots for each of magnitudes, overlaps, and stellar densities.

    :param hpx: HEALPIX value
    :param hpx_files: List of all HEALPIX values from present files
    """
    qrun = None

    # hpx_files is sorted, allowing a simple iteration over the list
    for filename in hpx_files:
        x = re.match('(?P<number>\d{3,5})_(?P<qrun>\d{2}[A-z]{2}\d{2})', filename)

        if qrun is None:
            qrun = x.group('qrun')

        if hpx in filename and qrun in filename and 'density' in filename:
            hdu = fits.open(FITS_DATA_DIR + filename)[0]
            w = wcs.WCS(hdu.header)
            plt.subplot(313, projection=w)
            plt.imshow(hdu.data, origin='lower', cmap='binary')
            plt.title(x.group('number') + ' ' + qrun + ' Stellar Densities')
            plt.xlabel('RA')
            plt.ylabel('Dec')

        elif hpx in filename and qrun in filename and '_mag_' in filename:
            hdu = fits.open(FITS_DATA_DIR + filename)[0]
            w = wcs.WCS(hdu.header)
            plt.subplot(311, projection=w)
            plt.imshow(hdu.data, origin='lower', cmap='viridis', vmin=np.min(hdu.data[np.nonzero(hdu.data)]) - 0.25)
            plt.title(x.group('number') + ' ' + qrun + ' Magnitudes')
            plt.xlabel('RA')
            plt.ylabel('Dec')

        elif hpx in filename and qrun in filename and 'overlap' in filename:
            hdu = fits.open(FITS_DATA_DIR + filename)[0]
            w = wcs.WCS(hdu.header)
            plt.subplot(312, projection=w)
            plt.imshow(hdu.data, origin='lower', cmap='binary')
            plt.title(x.group('number') + ' ' + qrun + ' Overlaps')
            plt.xlabel('RA')
            plt.ylabel('Dec')
            plt.show()  # overlap file comes last in the directory, wait to show the figure until it's reached
            qrun = None  # reset qrun to move onto next QRUNID set

        else:
            qrun = None


def load_stds(hpx, hpx_files):
    """
    Generates a histogram from the standard deviations of magnitudes in a given catalog (identified by its HEALPIX)

    :param hpx: HEALPIX value
    :param hpx_files: List of all HEALPIX values from present files
    """
    cat = storage.HPXCatalog(int(hpx), catalog_dir='catalogs/master', dest_directory='master')
    table = cat.table
    condition = (table['MATCHES'] > 2)
    hpxids = np.unique(table['HPXID'][condition])
    stds = None

    for filename in hpx_files:
        if str(hpx) in filename and 'stds' in filename:
            stds = fits.open(MAG_STD_DIR + filename)[1]
            break

    if stds is not None:
        plt.hist(stds.data.field('mag_std') / table['MAGERR_AUTO'][:len(hpxids)], bins=100, range=(-0.25, 10))
        # plt.hist([table['MAGERR_AUTO'][index] / i for index, i in enumerate(stds) if i != 0],
        #          50,
        #          range=(-0.2, 5))
        plt.xlabel("Magnitude St. Deviation / Magnitude Error ")
        plt.ylabel("Frequency")
        plt.title(str(hpx) + " Histogram of Magnitude Standard Deviation / MAGERR_AUTO")
        plt.show()
    else:
        print "Data files for HPX {} not written.".format(hpx)


def latitude_counts(hpx, hpx_files):

    for filename in sorted(hpx_files):
        if hpx in filename:
            x = re.match('(?P<hpx>\d{4}).*', filename)
            data = fits.open(HELIO_DATA_DIR + filename)[0].data
            plt.plot(data, '.')
            plt.xlabel("Latitude")
            plt.ylabel("Frequency")
            plt.title(x.group('hpx') + " Histogram of expected object counts")
            plt.show()


def main(params):
    if params.std_hist:
        direc = MAG_STD_DIR
        func = load_stds
    elif params.lat_count:
        direc = HELIO_DATA_DIR
        func = latitude_counts
    else:
        direc = FITS_DATA_DIR
        func = load_images

    try:
        directory = sorted(os.listdir(direc))
    except OSError as ex:
        raise ex

    hpx_values = []
    for name in directory:
        # HPX_02434_RA_185.6_DEC_+37.2
        x = re.match('(?P<number>\d{3,5})', name)
        if int(x.group('number')) not in hpx_values:
            hpx_values.append(int(x.group('number')))

    if len(hpx_values) == 0:
        print "Data files not found. Empty directory: {}".format(direc)
        sys.exit(0)

    while True:
        print "Enter one of: ",
        for i in sorted(hpx_values):
            print i,
        print
        print "Hit <enter> key to view all, or type 'exit' to quit."
        hpx = raw_input("Enter HPX to display: ")

        hpx_files = []
        for filename in directory:
            if hpx in filename:
                hpx_files.append(filename)

        if hpx == 'exit':
            sys.exit(0)

        elif hpx == '':
            for i in sorted(hpx_values):
                func(str(i), hpx_files)

        else:
            func(hpx, hpx_files)
