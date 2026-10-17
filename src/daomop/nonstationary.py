import argparse
import logging
from matplotlib import pyplot as plt
from numpy import *
from daomop import storage
from sklearn.ensemble import IsolationForest as IF
from scipy.spatial import cKDTree as KDTree
from daomop import util
import sys
import os
import pickle
DATADIR = __PATH__ = os.path.join(os.path.dirname(__file__), 'data')

class fits_cat(object):
    def __init__(self, hpx_catalog):
        self.catalog = hpx_catalog.table

        # Remove rows with infinite values from the table.
        self.array = self.feature(self.catalog)
        self.catalog = self.catalog[sum(isfinite(self.array), axis=1) == self.array.shape[1]]

        self.stationary0 = self.catalog[self.catalog['MATCHES'] > 0]
        self.stationary1 = self.catalog[self.catalog['MATCHES'] > 1]
        self.stationary2 = self.catalog[self.catalog['MATCHES'] > 2]
        self.ns0 = self.catalog[logical_and(logical_and(self.catalog['MATCHES'] == 0, self.catalog['OVERLAPS'] > 0),
                                            self.catalog['MAG_PSF'] < 27)]
        self.ns1 = self.catalog[logical_and(self.catalog['MATCHES'] == 0, self.catalog['OVERLAPS'] > 1)]
        self.ns2 = self.catalog[logical_and(self.catalog['MATCHES'] == 0, self.catalog['OVERLAPS'] > 2)]
        self.real_feature = self.feature(self.stationary2)
        self.ns_feature = self.feature(self.ns0)

    def plot(self):
        # plt.xlim(165.5, 168.0)
        # plt.ylim(31.5, 32.0)
        # plt.xlim(165, 168.0)
        # plt.plot(self.stationary0['X_WORLD'], self.stationary0['Y_WORLD'], '.', label='MATCH >0')
        # plt.plot(self.stationary1['X_WORLD'], self.stationary1['Y_WORLD'], '.', label='MATCH >1')
        # plt.plot(self.ns['X_WORLD'], self.ns['Y_WORLD'], '.', label='NS, OVERLAPS >1')
        # plt.plot(self.ns1['X_WORLD'], self.ns1['Y_WORLD'], '.', label='NS, OVERLAPS >2')
        plt.plot(self.catalog['X_WORLD'], self.catalog['Y_WORLD'], marker='.', label='ALL')
        plt.legend(loc=2)

    def feature(self, cat):
        return array([cat['MAG_ISO'], cat['MAGERR_ISO'], cat['MAG_ISOCOR'], cat['MAGERR_ISOCOR'],
                      cat['MAG_APER'], cat['MAGERR_APER'], cat['MAG_AUTO'], cat['MAGERR_AUTO'],
                      cat['MAG_PETRO'], cat['MAGERR_PETRO'], cat['MAG_WIN'], cat['MAGERR_WIN'],
                      cat['SNR_WIN'], cat['MAG_SOMFIT'], cat['MAGERR_SOMFIT'], cat['KRON_RADIUS'],
                      cat['PETRO_RADIUS'], cat['BACKGROUND'], cat['ISOAREA_IMAGE'], cat['ISOAREAF_IMAGE'],
                      cat['A_IMAGE'], cat['B_IMAGE'], cat['THETA_IMAGE'], cat['ERRX2_IMAGE'],
                      cat['ERRY2_IMAGE'], cat['ERRXY_IMAGE'], cat['ERRCXX_IMAGE'], cat['ERRCYY_IMAGE'],
                      cat['ERRCXY_IMAGE'], cat['CXXWIN_IMAGE'], cat['CYYWIN_IMAGE'], cat['CXYWIN_IMAGE'],
                      cat['AWIN_IMAGE'], cat['BWIN_IMAGE'], cat['ERRX2WIN_IMAGE'], cat['ERRY2WIN_IMAGE'],
                      cat['ERRXYWIN_IMAGE'], cat['ERRCXXWIN_IMAGE'], cat['ERRCYYWIN_IMAGE'], cat['ERRCXYWIN_IMAGE'],
                      cat['MU_THRESHOLD'], cat['MU_MAX'], cat['ISO0'], cat['ISO1'], cat['ISO2'], cat['ISO3'],
                      cat['ISO4'], cat['ISO5'], cat['ISO6'], cat['ISO7'], cat['FWHM_IMAGE'], cat['ELONGATION'],
                      cat['ELLIPTICITY'], cat['POLAR_IMAGE'], cat['POLARWIN_IMAGE'], cat['FLUX_RADIUS'],
                      cat['FWHMPSF_IMAGE'], cat['MAG_PSF'], cat['MAGERR_PSF'], cat['CHI2_PSF'],
                      cat['ERRCXXPSF_IMAGE'], cat['ERRCYYPSF_IMAGE'], cat['ERRCXYPSF_IMAGE'],
                      cat['ERRAPSF_IMAGE'], cat['ERRBPSF_IMAGE'], cat['ERRAPSF_WORLD'], cat['ERRBPSF_WORLD'],
                      cat['ERRTHETAPSF_SKY']]).T


class IF_real_bogus(object):
    def __init__(self, feature):
        self.feature = feature

    # self.Class = []
    # [self.Class.append('real') for i in range(len(self.feature))]
    # self.trainning, self.test, self.Class_train, self.Class_test = train_test_split(self.real_feature,\
    #                        self.real_Class, test_size=0.4, random_state=0)
    # print self.trainning
    def train(self):
        self.IFmod = IF(n_estimators=160)
        self.IFmod.fit(self.feature)

    def validation(self):
        # score = self.IFmod.decision_function(self.feature)
        result = self.IFmod.predict(self.feature)
        self.normal = (result == 1)
        self.abnormal = (result == -1)
        print "total:{0}, normal:{1}, abnormal:{2}, normal/total:{3}".format(len(self.feature),
                                                                             self.normal.sum(), self.abnormal.sum(),
                                                                             self.normal.sum() / float(
                                                                                 len(self.feature)))


def build_model(hpx_catalog):
    """
    Build a model of stellar sources to be used to select only stellar sources from the input catlogs.
    :param hpx_catalog:  The healpix catalog to use as the source of stellar profiles.
    :type hpx_catalog: storage.HPXCatalog
    :return: Model to be used for star galaxy separation.
    :rtype: IF_real_bogus.IFmod
    """

    logging.info("Building stationary source model based on stationary sources in {}".format(hpx_catalog.filename))
    sources = fits_cat(hpx_catalog)
    real_bogus = IF_real_bogus(sources.real_feature)
    real_bogus.train()
    real_bogus.validation()
    return real_bogus.IFmod


def remove_ps1_source(ps1_cat, ns_cat):
    ps1 = ps1_cat.table
    ps1_ra = ps1['f0']
    ps1_dec = ps1['f1']
    ra = ns_cat['X_WORLD']
    dec = ns_cat['Y_WORLD']
    ns_tree = KDTree(zip(ra, dec))
    ps1_tree = KDTree(zip(ps1_ra, ps1_dec))
    match = ns_tree.query_ball_tree(ps1_tree, r=.5 / 3600.)
    # print len(match), len(ra), len(ps1)
    ps1_match = []
    for i in match:
        if i != []:
            ps1_match.append(True)
        else:
            ps1_match.append(False)
    return array(ps1_match)


def ns_cat(hpx_catalog, ps1_cat, rb_model):

    sources = fits_cat(hpx_catalog)
    real_bogus = IF_real_bogus(sources.ns_feature)
    real_bogus.IFmod = rb_model
    real_bogus.validation()
    ns_cat = sources.ns0[real_bogus.normal]
    ps1_match = remove_ps1_source(ps1_cat, ns_cat)

    logging.info("PS1 remove {0} stationary sources".format(ps1_match.sum()))
    logging.info("{0} non-stationary sources remain".format((~ps1_match).sum()))
    nonstationary = storage.HPXCatalog(hpx_catalog.pixel,
                                       version='_ns_cat',
                                       catalog_dir=hpx_catalog.catalog_dir)
    nonstationary.hdulist.append(hpx_catalog.hdulist[0])
    logging.info("Writing non-stationary catalog to {}".format(nonstationary.filename))
    nonstationary.table = ns_cat[~ps1_match]
    nonstationary.write()
    nonstationary.put()


def main():
    """
    Retrieve HPX based catalogs from VOSpace storage, build a stationary source model from those and then
    remove the stationary sources from the input catalogs.
    :param training_hpx:
    :param search_hpx:
    :return:
    """

    parser = argparse.ArgumentParser(
        description='Creates a catalog of non-stationary sources based on input catalogs and NN.')
    parser.add_argument("training_hpx", type=int,
                        help="Healpix to use as the training set input.")
    parser.add_argument("--ps1_dir",
                        default="ps_pv3_STK_hp",
                        help="Catalog directory that contains the PS1 catalogs")
    parser.add_argument("--dbimages",
                        action="store",
                        default="vos:cfis/solar_system/dbimages",
                        help='vospace dbimages containerNode')
    parser.add_argument("--catalogs",
                        action="store",
                        default="catalogs",
                        help='dbimages subdirectory where catalogs will be stored.')
    parser.add_argument("healpix",
                        type=int,
                        help="healpix to process")
    parser.add_argument("--dry-run",
                        action="store_true",
                        help="DRY RUN, don't copy results to VOSpace, implies --force")
    parser.add_argument("--verbose", "-v",
                        action="store_true")
    parser.add_argument("--force", default=False,
                        action="store_true")
    parser.add_argument("--debug", "-d",
                        action="store_true")
    parser.add_argument("qrunid", help="The CFHT QRUN to build non-stationary catalogs for.")

    cmd_line = " ".join(sys.argv)
    args = parser.parse_args()
    util.set_logger(args)
    logging.info("Started {}".format(cmd_line))

    storage.DBIMAGES = args.dbimages

    catalog_dirname = "{}/{}".format(args.catalogs, args.qrunid)

    training_catalog = storage.HPXCatalog(pixel=args.training_hpx, catalog_dir=catalog_dirname)
    stellar_model = storage.HPXCatalog(pixel=args.training_hpx, catalog_dir=catalog_dirname,
                                       ext=".pkl", version="_IFmod")
    # Get or build a model of the stellar profile.
    try:
        stellar_model.get()
        with open(stellar_model.filename, 'r') as fobj:
            rb_model = pickle.load(fobj)
    except Exception as ex:
        logging.warning(str(ex))
        rb_model = build_model(training_catalog)
        with open(stellar_model.filename, 'w') as fobj:
            pickle.dump(rb_model, fobj)
        stellar_model.put()

    ps1_catalog = storage.HPXCatalog(pixel=args.healpix, catalog_dir=args.ps1_dir, version="_STK")
    search_catalog = storage.HPXCatalog(pixel=args.healpix, catalog_dir=catalog_dirname)
    ns_cat(search_catalog, ps1_catalog, rb_model)


if __name__ == '__main__':
    main()
