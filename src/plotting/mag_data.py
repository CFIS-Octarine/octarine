"""
FITS data file generator. Writes sky coverage plots from source master catalogs on VOSpace. Generates 3 unique files
for every QRUNID in each master catalog (unless the data files are empty, in which only a single file is written to
signify zero data).

Defaulted to generate data files from every master catalog.
In main(), alter the line
    > if hpx not in reg and hpx > 0:
to be
    > if hpx not in reg and hpx == <HEALPIX>:
to generate files for a single catalog.

Usage:
~/octarine/src$ python plot.py --mag
"""

import os
import numpy as np
import re
from astropy.io import fits
from astropy import wcs

from daomop import storage

PIXEL_WIDTH = 90.0/3600.0
FITS_DIRECTORY = 'plotting/fits_data/'


class IDX(object):
    def __init__(self, range_min, range_max, width):
        self.min = range_min
        self.max = range_max
        self.width = width
        self.nbin = int((range_max - range_min) / width) + 1
        self.centre = (self.max+self.min) / 2

    def __len__(self):
        return self.nbin


def fits_factory(hpx):
    """
    Generate FITS data files from master catalog on VOSpace. Three files will be created per QRUNID per catalog.
    Each file is created via binning the XWORLD/YWORLD sky plots, with the condition that the flux limit > 0.15.
    Each bin is 90"x90".

    The first data file will be a coverage plot where each pixel contains the third faintest magnitude for that bin (
    this requires >2 data points in the pixel).
    The second file will contain the amount of object detected in the bin.
    The third file will be the stellar density, and does not require there to be >2 data points in the bin.

    In the case where there are no satisfying data points for the first two files, only the stellar density file will
    be written.

    :param hpx: HEALPIX value to specify which catalog to build images from
    """
    try:
        print hpx

        cat = storage.HPXCatalog(hpx, catalog_dir='catalogs/master', dest_directory='master')
        table = cat.table

        table['frame'] = [dataset_name.split('p')[0] for dataset_name in table['dataset_name']]
        dataset_names = np.unique(table['frame'])
        table_entries = [table[table['frame'] == name] for name in dataset_names]

        mag_lims = {}

        for table_entry in table_entries:
            previous_mag = None
            for mag in np.arange(table_entry['MAG_AUTO'].min(), table_entry['MAG_AUTO'].max(), 0.20):
                condition = np.all((table_entry['MAG_AUTO'] > mag, table_entry['MAG_AUTO'] <= mag + 0.2), axis=0)
                magerr = np.median(table_entry['MAGERR_AUTO'][condition])
                if magerr > 0.15:
                    if previous_mag is None:
                        previous_mag = mag
                    else:
                        frame = table_entry['frame'][0].split('p')[0]
                        mag_lims[frame] = mag
                        break
                else:
                    previous_mag = None

        ra_idx = IDX(table['X_WORLD'].min(), table['X_WORLD'].max(), PIXEL_WIDTH)
        dec_idx = IDX(table['Y_WORLD'].min(), table['Y_WORLD'].max(), PIXEL_WIDTH)

        mag_data = np.zeros((len(dec_idx), len(ra_idx)))

        w = wcs.WCS(naxis=2)
        w.wcs.cd = [[PIXEL_WIDTH, 0], [0, PIXEL_WIDTH]]
        crpix1 = mag_data.shape[1]/2
        crpix2 = mag_data.shape[0]/2
        w.wcs.crpix = [crpix1, crpix2]
        crval = [ra_idx.centre, dec_idx.centre]
        w.wcs.crval = crval
        w.wcs.ctype = ['ra--tan', 'dec-tan']
        w.wcs.cunit = ['deg', 'deg']
        header = w.to_header()

        print mag_data.shape[1], mag_data.shape[0]
        for qrun in np.unique(table['QRUNID']):

            # set file names and reset data arrays for each qrunid
            if len(str(hpx)) == 3:
                hpx = '0' + str(hpx)
            mag_image_filename = FITS_DIRECTORY + str(hpx) + '_' + qrun + '_mag_data.fits'
            overlap_image_filename = FITS_DIRECTORY + str(hpx) + '_' + qrun + '_overlap_image.fits'
            density_image_filename = FITS_DIRECTORY + str(hpx) + '_' + qrun + '_density_image.fits'
            mag_data = np.zeros((len(dec_idx), len(ra_idx)))
            overlap_data = np.zeros((len(dec_idx), len(ra_idx)))
            stellar_density_data = np.zeros((len(dec_idx), len(ra_idx)))
            print qrun, mag_image_filename

            # if the density file exists, this set has been examined/written before and can be skipped
            if not os.path.exists(density_image_filename):
                count = 0
                for xx in range(mag_data.shape[1]):
                    print count
                    count += 1

                    for yy in range(mag_data.shape[0]):
                        ra, dec = w.all_pix2world(xx, yy, 0)

                        ra_cond = np.all((table['X_WORLD'] >= ra,
                                          table['X_WORLD'] < ra + PIXEL_WIDTH,
                                          table['OVERLAPS'] > 1,
                                          table['QRUNID'] == qrun),
                                         axis=0)

                        dec_cond = np.all((ra_cond,
                                           table['Y_WORLD'] >= dec,
                                           table['Y_WORLD'] < dec + PIXEL_WIDTH),
                                          axis=0)

                        dataset_names = np.unique(table[dec_cond]['frame'])

                        if len(dataset_names) != 0:
                            density = float(len(table[dec_cond])) / len(dataset_names)
                            stellar_density_data[yy, xx] = density

                            if len(dataset_names) > 2:
                                l = []
                                try:
                                    for frame in dataset_names:
                                        l.append(float(mag_lims[frame]))  # getting third faintest magnitude
                                    third_faintest = sorted(l)[-3]
                                    print 'index: ({} , {})'.format(yy, xx)
                                    mag_data[yy, xx] = third_faintest
                                    overlap_data[yy, xx] = len(dataset_names)
                                except KeyError:
                                    print KeyError
                                    continue

                for row in mag_data:
                    if row.any() != 0:
                        # only write these files if there's data that has been inserted
                        mag_image = fits.PrimaryHDU(data=mag_data, header=header)
                        mag_image.writeto(mag_image_filename)

                        overlap_image = fits.PrimaryHDU(data=overlap_data, header=header)
                        overlap_image.writeto(overlap_image_filename)
                        break

                # write the stellar density image regardless so it's not regathered again in case of empty mag_data
                density_image = fits.PrimaryHDU(data=stellar_density_data, header=header)
                density_image.writeto(density_image_filename)

    except Exception as ex:
        print ex
        return


def main():
    if not os.path.exists(FITS_DIRECTORY):
        os.mkdir(FITS_DIRECTORY)

    directory = storage.listdir(os.path.join(os.path.dirname(storage.DBIMAGES),
                                             storage.CATALOG, 'master'), force=True)
    reg = []
    for item in directory:
        # HPX_02434_RA_185.6_DEC_+37.2
        x = re.match('HPX_(?P<pix>\d{5})_RA_(?P<ra>\d{3}\.\d)_DEC_\+(?P<dec>\d{2}\.\d)', item)

        hpx = int(x.group('pix'))

        if hpx not in reg and hpx > 0:  # alter 'hpx > 0' to use specific files
            reg.append(hpx)
            fits_factory(hpx)
