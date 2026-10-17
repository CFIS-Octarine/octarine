"""Build the source detection list from an image.

This source detection system uses sextactor and produces PSF based measurement values."""
import argparse
import logging
import os
import subprocess
import sys
import traceback
from . import storage
from . import util

task = 'build_cat'

__PATH__ = os.path.dirname(__file__)
SEX_CONFIG = os.path.join(__PATH__, 'config')
os.environ['SEX_CONFIG'] = SEX_CONFIG


def run(expnum, ccd, version, prefix, dry_run, force):

    message = storage.SUCCESS

    observation = storage.Observation(expnum)
    image = storage.FitsImage(observation, ccd=ccd)
    if image.status(task) and not force:
        logging.info("{} completed successfully for {} {} {} {}".format(task, prefix, expnum, version, ccd))
        return

    with storage.LoggingManager(task, prefix, expnum, ccd, version, dry_run):
        try:
            image.get(return_file=True, convert_to_sip=False)
            image.flat_field.get(return_file=True, convert_to_sip=False)

            # Build the PSF model input catalog
            logging.info("Building PSF input catalog")
            logging.info("Using config: {}".format(os.path.join(SEX_CONFIG, 'pre_psfex.sex')))
            ldac_catalog = storage.Artifact(observation, ccd=ccd, ext=".ldac")
            sex_cmd = 'sex'
            for filename in ['/opt/local/bin/sex', '/usr/bin/sex', '/usr/local/bin/sex']:
                if os.access(filename, os.X_OK):
                    sex_cmd = filename
                    break
            cmd = [sex_cmd, image.filename,
                   '-c', os.path.join(SEX_CONFIG, 'pre_psfex.sex'),
                   '-CATALOG_NAME', ldac_catalog.filename,
                   '-WEIGHT_IMAGE', image.flat_field.filename,
                   '-MAG_ZEROPOINT', str(image.zeropoint)]
            logging.info(" ".join(cmd))
            logging.info(subprocess.check_output(cmd, stderr=subprocess.STDOUT))

            # Build the PSF model
            cmd = ['psfex', ldac_catalog.filename,
                   '-c', os.path.join(SEX_CONFIG, 'default.psfex')]
            logging.info(" ".join(cmd))
            logging.info(subprocess.check_output(cmd,
                                                 stderr=subprocess.STDOUT))

            # Build a source catalog using the PSF model.
            fits_catalog = storage.Artifact(observation, ccd=ccd, ext=".cat.fits")
            psf = storage.Artifact(observation, ccd=ccd, ext=".psf")
            cmd = [sex_cmd,
                   '-c', os.path.join(SEX_CONFIG, 'ml.sex'),
                   '-WEIGHT_IMAGE', image.flat_field.filename,
                   '-CATALOG_NAME', fits_catalog.filename,
                   '-PSF_NAME', psf.filename,
                   '-MAG_ZEROPOINT', str(image.zeropoint),
                   image.filename]
            logging.info(" ".join(cmd))
            logging.info(subprocess.check_output(cmd, stderr=subprocess.STDOUT))

            if dry_run:
                return

            # transfer results to storage.
            fits_catalog.put()
            psf.put()
            logging.info(message)

        except Exception as e:
            logging.debug(traceback.format_exc())
            logging.error(str(e))
            message = str(e)

    image.status(task, message)

    return


def main():

    parser = argparse.ArgumentParser(
        description='Run biuld_cat chunk of the CFIS pipeline')

    parser.add_argument('--ccd', '-c',
                        action='store',
                        type=int,
                        dest='ccd',
                        default=None,
                        help='which ccd to process, default is all')
    parser.add_argument("--dbimages",
                        action="store",
                        default="vos:cfis/solar_system/dbimages",
                        help='vospace dbimages containerNode')
    parser.add_argument("expnum",
                        type=int,
                        nargs='+',
                        help="expnum(s) to process")
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

    prefix = ''
    version = 'p'

    storage.DBIMAGES = args.dbimages

    exit_code = 0
    for expnum in args.expnum:
        if args.ccd is None:
            if int(expnum) < 1785619:
                # Last exposures with 36 CCD MegaPrime
                ccdlist = range(0, 36)
            else:
                # First exposures with 40 CCD MegaPrime
                ccdlist = range(0, 40)
        else:
            ccdlist = [args.ccd]
        for ccd in ccdlist:
            run(expnum, ccd, version, prefix, args.dry_run, args.force)
    return exit_code

if __name__ == '__main__':
    sys.exit(main())
