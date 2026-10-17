import argparse
import sys
import logging
import errno
from . import util
from . import storage
from .storage import archive_url
from .storage import Observation
from .storage import FitsImage, Header
from .storage import make_path
from .storage import make_link
from .storage import pitcairn_uri
from .storage import isfile


def run(dataset_name):
    """Given a dataset_name created the desired dbimages directories
    and links to the raw and processed data files stored at CADC and vospace.

    @param dataset_name: the target_name of the CFHT dataset to make a link to.
    """

    observation = Observation(dataset_name)

    version = 'o'
    ext = '.fits.fz'
    artifact = FitsImage(observation, version=version, ext=ext)
    source = archive_url(dataset_name, version)
    logging.debug("Making link between {} and {}".format(source, artifact.uri))
    make_path(artifact.uri)
    make_link(source, artifact.uri)

    logging.debug("Linking to RAW header")
    source = archive_url(dataset_name, version=version, fhead='true')
    header = Header(observation, version=version)
    make_link(source, header.uri)

    logging.debug("Making link between PROC'd image and dbimages")
    version = 'p'
    ext = '.fits.fz'
    artifact = FitsImage(observation, version=version, ext=ext)
    # Source is either PITCAIRN processing or CFHT archive URL
    source = pitcairn_uri(dataset_name)
    if not isfile(source):
        # Try pitcairn without the .fz
        source = pitcairn_uri(dataset_name, ext=".fits")
    if not isfile(source):
        raise OSError(errno.EEXIST, "No file in pitcairn to link to for {}".format(source), dataset_name)

    make_link(source, artifact.uri)

    logging.debug("Link up the header")
    # Can be the CFHTSG header or the one in the CFHT archive one.
    source = archive_url(dataset_name, version, ext='.head', archive='CFHTSG')
    header = Header(observation, version=version)
    if not isfile(source):
        raise OSError(errno.EEXIST, "No header to link to for {}".format(source), dataset_name)

    make_link(source, header.uri)

    return True


def main():
    parser = argparse.ArgumentParser(
        description='Create the vospace entries required for pipeline processing')

    parser.add_argument("--dbimages",
                        action="store",
                        default="vos:cfis/solar_system/dbimages",
                        help='vospace dbimages containerNode')
    parser.add_argument("expnum",
                        type=int,
                        nargs='+',
                        help="expnum(s) to create directories for")
    parser.add_argument("--dry-run",
                        action="store_true",
                        help="DRY RUN, don't copy results to VOSpace, implies --force")
    parser.add_argument("--verbose", "-v",
                        action="store_true")
    parser.add_argument("--force", default=False,
                        action="store_true")
    parser.add_argument("--debug", "-d",
                        action="store_true")
    parser.add_argument("--pitcairn", default="vos:cfis/pitcairn", action="store", help="vospace containing pitcairn processed images")

    cmd_line = " ".join(sys.argv)
    args = parser.parse_args()

    util.set_logger(args)
    logging.info("Started {}".format(cmd_line))

    storage.DBIMAGES = args.dbimages
    storage.PITCAIRN = args.pitcairn

    exit_code = 0
    for expnum in args.expnum:
        run(expnum)
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
