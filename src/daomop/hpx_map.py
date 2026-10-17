import matplotlib as mpl
mpl.use('Agg')
import healpy as hp
import numpy as np
from astropy.table import Table
from matplotlib import pyplot
import os
import sys
import re
import logging
import argparse

from . import storage
from . import util 

def run(healpix, min_matches=0, max_matches=10, min_overlaps=0, max_overlaps=10, dry_run=False):
    """
    Plot the RA/DEC values of all objects in a given file, colour points by value of expnum.
    
    :param: healpix to build coverage map for
    """

    logging.debug("Building coverage map using healpix: {}".format(healpix))
    catalog = storage.HPXCatalog(healpix)
    t = catalog.table

    expnums = np.array([x.split('p')[0] for x in t['dataset_name']])
    colours = ['r', 'g', 'b', 'y']
    count = 0

    for expnum in np.unique(expnums):
        colour = colours[count % len(colours)]
        cond = np.all((t['FLUX_RADIUS'] > 2.5, expnums == expnum, t['MATCHES'] >= min_matches, t['MATCHES'] <= max_matches, t['OVERLAPS'] >= min_overlaps, t['OVERLAPS'] <= max_overlaps), axis=0)
        pyplot.plot(t['X_WORLD'][cond], t['Y_WORLD'][cond], ',{}'.format(colour), ms=1, alpha=.25, label=str(expnum))
        count += 1

    pyplot.legend(fontsize=8)
    pyplot.xlabel("RA")
    pyplot.ylabel("DEC")
    pyplot.title(catalog.filename)

    # Now convert to phi,theta representation:
    phi_theta = util.healpix_to_corners(healpix)
    xy = list(phi_theta)
    xy.append(phi_theta[0])
    xy = np.array(xy).transpose()
    pyplot.plot(xy[0], xy[1], '-k')
    artifact = storage.HPXCatalog(healpix, ext=".png")
    pyplot.savefig(artifact.filename)
    if not dry_run:
        artifact.put()


def main():
    parser = argparse.ArgumentParser(
        description='Create a sky plot of the sources detected in the given HEALPIX.')

    parser.add_argument("--dbimages",
                        action="store",
                        default="vos:cfis/solar_system/dbimages",
                        help='vospace dbimages containerNode')
    parser.add_argument("healpix",
                        type=int,
                        help="healpix to map")
    parser.add_argument("--min-matches",
                        type=int,
                        action="store",
                        default=0,
                        help="Minimum number of matches a source must have to be included in the map")
    parser.add_argument("--min-overlaps",
                        type=int,
                        action="store",
                        default=0,
                        help="Minimum number of overlaping images a source must have to be included in the map")
    parser.add_argument("--max-matches",
                        default=100,
                        action="store",
                        type=int,
                        help="Maximum number of matches a source must have to be included in the map")
    parser.add_argument("--max-overlaps",
                        default=100,
                        action="store",
                        type=int,
                        help="Maximum number of overlaping images a source can have and be include in the map")
    parser.add_argument("--dry-run",
                        action="store_true",
                        help="DRY RUN, don't copy results to VOSpace, implies --force")
    parser.add_argument("--verbose", "-v",
                        action="store_true")
    parser.add_argument("--force", default=False,
                        action="store_true")
    parser.add_argument("--debug", "-d",
                        action="store_true")

    cmd_line = " ".join(sys.argv)
    args = parser.parse_args()

    util.set_logger(args)
    logging.info("Started {}".format(cmd_line))

    storage.DBIMAGES = args.dbimages
    prefix = ''
    version = 'p'

    exit_code = 0
    run(args.healpix, min_matches=args.min_matches, max_matches=args.max_matches, min_overlaps=args.min_overlaps, max_overlaps=args.max_overlaps, dry_run=args.dry_run)
    return exit_code

if __name__ == '__main__':
    sys.exit(main())
