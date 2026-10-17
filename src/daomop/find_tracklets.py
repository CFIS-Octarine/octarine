import logging
import os
import sys
import time
import copy
import json
import warnings
import argparse
import numpy
from astropy.time import Time
from astropy import units
from astropy.coordinates import SkyCoord
import numpy as np
from matplotlib import pyplot

from . import storage
from . import util

warnings.filterwarnings("ignore")

__version__ = "2.4"


deduplicate = True
# the setup for multiprocess
multi = False
N_threading = 1

# the paramaters for BK diagnosis
printmpc = "yes"
obscode = "568"
crit_residual = 0.3
maxau = 100.
minau = 30.
n_proc = 3
astrometryerror = 1. / 3600  # arcsec

# ProcessFilterid = ['HSC-g', 'HSC-r', 'HSC-i', 'HSC-r2', 'HSC-i2']
Processmag = 'MAG_PSF'
Processmagerr = 'MAGERR_PSF'
ProcessMinExptime = 10.
stackfitsname = 'allmlns.fits'
allpairname = 'allpairs.json'
alltrackname = 'mjdalltracks.json'


#####################################
# Searching moving objects
#####################################
class SearchMovingObjects:

    def __init__(self, non_stationary_catalog):
        """

        :param non_stationary_catalog: HPX calatog from DAOMOP with source detections.
        :type non_stationary_catalog: storage.HPXCatalog
        :param rates:
        """
        self.catalog = non_stationary_catalog
        self.catalog.table['DATE'] = Time(self.catalog.table['mid_mjdate'], format='mjd')
        self.catalog.table['SKYCOORD'] = SkyCoord(self.catalog.table['X_WORLD'],
                                                 self.catalog.table['Y_WORLD'],
                                                 unit='degree')

        self.midmjd = numpy.median(numpy.unique(self.catalog.table['mid_mjdate']))
        self.coordinate = storage.util.healpix_to_skycoord(self.catalog.pixel)
        self.cra = self.coordinate.ra.degree
        self.cdec = self.coordinate.dec.degree
        self.rate_low = 0.5 * units.arcsec/units.hour
        self.rate_high = 15 * units.arcsec/units.hour
        self.merr_limit = 0.25

        # the field/healpix radius, in degree [2.5 for moving rate up to 75 degree, 2 for 60 degree, 1.3 for < 45]
        self.fr = 2.5
        # the search step for searchline, in arcsecond
        self.step = 2.0
        # the search radius for searchline, in arcsecond
        self.sr = 5.0

        # Width of the search cone
        self.openangle = 30.
        # Angle of search cone
        self.avgang = -90
        # step size in search space
        self.angle_step = 1
        # the fastet moving speed for the object at * au, the tracks faster than this value will be ignored.
        fastau = 3.0
        slowau = 30.0  # the slow moving speed for the object at * au, slower tracks *LOOSen CONDITION OF SPEED *
        N_dets = 3  # only consider the tracks with N_dets at least
        max_dm = 1.0  # only consider the tracks with dm < max_dm
        min_motion = 0.2  # in arcsecond, the requirement of minmum motion between two exposures
        min_dt = 20. / 60 / 24  # in day, the requirement of minmum time separation.
        max_v = 15.0 * 24 / 3600  # degree per day
        min_v = 1.0 * 24 / 3600  # degree per day

    def searchline(self, inputra, inputdec, ang, sr):
        """
        Find all values in catalog where the RA/DEC coordinate
        :param inputra:
        :param inputdec:
        :param ang:
        :param sr:
        :return:
        """

        slope = np.tan(np.radians(ang))
        b = inputdec - slope * inputra
        dr = (self.catalog.table['X_WORLD'] - inputra) - slope*(self.catalog.table['Y_WORLD'] - inputdec)
        snsra, snsdec = self.sns['X_WORLD'], self.sns['Y_WORLD']
        inline = abs(slope * snsra - snsdec + b) / ((slope ** 2 + 1) ** 0.5) < sr / 3600.0

        if sum(inline) == 0:
            return 0
        else:
            match = self.sns[inline]
            return match

    def oneangle(self, ang):
        t2 = time.time()
        result2 = []
        logging.info('[%s] Processing angle: %s' % (time.strftime("%D %H:%M:%S"), ang))
        for sl, ldec in enumerate(np.arange(self.cdec - self.fr, self.cdec + self.fr, self.step / 3600)):
            inputra, inputdec = self.cra, ldec
            mat = self.searchline(inputra, inputdec, ang, sr)
            if mat != 0:
                pairs = self.findpair(sl, mat, ang, fastau, slowau)
                result2.append(pairs)
        print '[%s] Processing angle: %s. Total pair: %s. Total time: %s' % (
        time.strftime("%D %H:%M:%S"), ang, len(result2), time.time() - t2)
        return result2

    def findpair(self, sl, mat, ang, fastau, slowau):
        vmax = 147. / fastau * 24 / 3600
        vmaxabs = abs(vmax)
        vslow = 147. / slowau * 24 / 3600.
        vslowabs = abs(vslow)
        line = mat
        # mjdlist = list(set([int(i - 2400000.5) for i in line['jd']]))
        mjdlist = list(set([int(i) for i in line['mid_mjdate']]))
        pairs = {}
        # 		for m in mjdlist:
        # 			#mdata = line[abs(line['mid_mjdate'] - m) < 1.0]
        # 			mdata = line[(line['mid_mjdate'] - m < 1.0)&(line['mid_mjdate'] - m > 0)]
        # 			mjdmax, mjdmin = mdata['mid_mjdate'].max(), mdata['mid_mjdate'].min()
        if True:
            m = mjdlist[0]
            mdata = line
            linking = {}
            for n, i in enumerate(sorted(mdata, key=lambda mdata: np.degrees(mdata['X_WORLD']))):
                tra, tdec, tmjd = float(i['X_WORLD']), float(i['Y_WORLD']), float(i['mid_mjdate'])
                tmag = float(i[Processmag])
                tmagerr = float(i[Processmagerr])
                tA = float(i['A_IMAGE'])
                tB = float(i['B_IMAGE'])
                ttheta = float(i['THETA_IMAGE'])
                tpath, tf, tet = i['dataset_name'], 'r', 30.0
                radecstr = '%13.9f_%13.11f_%03i_%04i_%04i' % (tra, tdec, ang, sl, n)
                linking[radecstr] = {'ra': [tra], 'dec': [tdec], 'mjd': [tmjd], 'mag': [tmag], 'magerr': [tmagerr],
                                     'A': [tA], 'B': [tB], 'theta': [ttheta], 'fitsname': [tpath], 'filterid': [tf],
                                     'exptime': [tet]}
                for obj in linking.keys():
                    tk = linking[obj]
                    if tmjd != tk['mjd'][0] and tmjd != tk['mjd'][-1]:
                        a_dt = abs(tmjd - tk['mjd'][-1])
                        d2 = ((tdec - tk['dec'][0]) ** 2 + (tra - tk['ra'][0]) ** 2) ** 0.5
                        a_v = d2 / a_dt
                        if a_dt > min_dt and a_v < max_v and a_v > min_v:
                            motion = ((tdec - tk['dec'][-1]) ** 2 + (tra - tk['ra'][-1]) ** 2) ** 0.5
                            if motion > min_motion / 3600.:
                                # going to add detection after detection, no matter prograde or retrograde
                                v2 = d2 / (tmjd - tk['mjd'][0])
                                v2abs = abs(v2)
                                newobj = '%s_%04i' % (obj, n)
                                if abs(v2) < abs(vmax) and abs(tk['mag'][0] - tmag) < max_dm:
                                    if tk['ra'].__len__() == 1:
                                        if tk['ra'][0] != tra:
                                            linking[newobj] = copy.deepcopy(tk)
                                            linking[newobj]['ra'].append(tra)
                                            linking[newobj]['dec'].append(tdec)
                                            linking[newobj]['mjd'].append(tmjd)
                                            linking[newobj]['mag'].append(tmag)
                                            linking[newobj]['magerr'].append(tmagerr)
                                            linking[newobj]['A'].append(tA)
                                            linking[newobj]['B'].append(tB)
                                            linking[newobj]['theta'].append(ttheta)
                                            linking[newobj]['fitsname'].append(tpath)
                                            linking[newobj]['filterid'].append(tf)
                                            linking[newobj]['exptime'].append(tet)
                                    else:
                                        if (tk['mjd'][-1] - tk['mjd'][0]) * (tmjd - tk['mjd'][-1]) > 0:
                                            goodvelocity = False
                                            # The condition for the larger error ratio of slow
                                            # moving objects in short duration
                                            d1 = ((tk['dec'][-1] - tk['dec'][0]) ** 2 + (
                                            tk['ra'][-1] - tk['ra'][0]) ** 2) ** 0.5
                                            v1 = d1 / (tk['mjd'][-1] - tk['mjd'][0])
                                            v1abs = abs(v1)
                                            # condition for slow mover (intra-night)
                                            if v1abs < vslowabs and v2abs < vslowabs:
                                                # if abs( (v2-v1)/min([v1, v2]) ) < 1 and abs(v1) < vmaxabs:
                                                goodvelocity = True
                                            # condition for fast mover (intra-night)
                                            else:
                                                if abs((v2 - v1) / min([v1abs, v2abs])) < 0.2 and v1abs < vmaxabs:
                                                    goodvelocity = True
                                            if goodvelocity:
                                                linking[newobj] = copy.deepcopy(tk)
                                                linking[newobj]['ra'].append(tra)
                                                linking[newobj]['dec'].append(tdec)
                                                linking[newobj]['mjd'].append(tmjd)
                                                linking[newobj]['mag'].append(tmag)
                                                linking[newobj]['magerr'].append(tmagerr)
                                                linking[newobj]['A'].append(tA)
                                                linking[newobj]['B'].append(tB)
                                                linking[newobj]['theta'].append(ttheta)
                                                linking[newobj]['fitsname'].append(tpath)
                                                linking[newobj]['filterid'].append(tf)
                                                linking[newobj]['exptime'].append(tet)

            # To remove the duplicate tracks and subtracks from matched pairs.
            # linkingout = copy.deepcopy(linking)
            if deduplicate:
                dlist = []
                for obj1 in linking.keys():
                    if linking[obj1]['ra'].__len__() < N_dets or max(linking[obj1]['mag']) - min(
                            linking[obj1]['mag']) > max_dm:
                        # del linkingout[obj1]
                        dlist.append(obj1)
                for obj1 in dlist:
                    del linking[obj1]
                # linking = copy.deepcopy(linkingout)
                dlist = []
                for obj1 in linking.keys():
                    for obj2 in linking.keys():
                        if obj1 != obj2:
                            if set(linking[obj1]['ra']) < set(linking[obj2]['ra']):
                                # del linkingout[obj1]
                                dlist.append(obj1)
                                break
                for obj1 in dlist:
                    del linking[obj1]
            # pairs.append(linkingout)
            # pairs[m] = linkingout
            pairs[m] = linking
        return pairs

    def search(self):

        for idx1 in range(len(self.catalog.table)):
            first_source = self.catalog.table[idx1]
            rate = first_source['SKYCOORD'].separation(self.catalog.table['SKYCOORD']) / (
                first_source['DATE'] - self.catalog.table['DATE'])
            # select the sources whose location is compatible with being linked.
            selected = np.all((rate < -self.rate_low,
                               rate > -self.rate_high,
                               self.catalog.table['MAGERR_AUTO'] < self.merr_limit,
                               self.catalog.table['DATE'] >= first_source['DATE']), axis=0)
            logging.debug("{} sources within rate cuts of {}".format(np.sum(selected),
                                                                    self.catalog.table['SKYCOORD'][idx1]))
            if not len(np.unique(self.catalog.table['DATE'][selected])) > 1:
                continue
            print self.catalog.table['SKYCOORD', 'DATE'][selected]
            rate1 = rate[selected]
            rate1.shape = (1, len(rate1))
            rate_matrix = np.abs((rate1 - rate1.T)/rate1)
            print rate_matrix
            rate_matrix = rate_matrix < 0.1
            row_count = 0
            print rate_matrix
            for row in rate_matrix:
                row_count += 1
                in_bounds = row * (np.arange(len(row)) > row_count)
                print in_bounds
                if not sum(in_bounds) > 0:
                    continue
                these_sources = self.catalog.table[selected][in_bounds]
                second_source = these_sources[0]
                print these_sources['SKYCOORD', 'DATE']
                for third_source in these_sources[these_sources['DATE'] > second_source['DATE']]:
                    a = 1.0*units.hour * (second_source['SKYCOORD'].separation(third_source['SKYCOORD'])) / (
                        second_source['DATE'] - third_source['DATE'])
                    b = 1.0*units.hour * (first_source['SKYCOORD'].separation(second_source['SKYCOORD'])) / (
                        first_source['DATE'] - second_source['DATE'])
                    c = 1.0*units.hour * (first_source['SKYCOORD'].separation(third_source['SKYCOORD'])) / (
                        first_source['DATE'] - third_source['DATE'])
                    A = np.arccos(np.cos(a) - np.sin(c)*np.sin(b))/(np.sin(c)*np.sin(b))

                    logging.info("Rates and opening angle: {} {} {} {}".format(a, b, c, A))
                    pyplot.plot(first_source['X_WORLD'], first_source['Y_WORLD'], 'or')
                    pyplot.plot(second_source['X_WORLD'], second_source['Y_WORLD'], 'og')
                    pyplot.plot(third_source['X_WORLD'], third_source['Y_WORLD'], 'ob')
                    pyplot.text((third_source['X_WORLD'] + first_source['X_WORLD'])/2.0,
                                (third_source['Y_WORLD'] + first_source['Y_WORLD']) / 2.0,
                                numpy.degrees(A))

                    pyplot.show()
        return


    def CleanTracks_array(self):
        mjdalltracks = {}
        tracksarray = np.array([], dtype=[('ra', 'f8'), ('dec', 'f8'), ('mjd', 'i4'), ('objid', 'S150')])
        for n, allsearchang in enumerate(s.results):
            t2 = time.time()
            angtrackslist = s.results[n]
            for i, alltrack in enumerate(angtrackslist):
                angtracks = angtrackslist[i]
                for mjd in angtracks.keys():
                    for obj in angtracks[mjd]:
                        # Checking if mjd is in mjdalltracks dictionary
                        if mjd in mjdalltracks.keys():
                            # Checking the deplicates
                            ifsub = []
                            verify = tracksarray[tracksarray['mjd'] == int(mjd)]
                            verify = verify[abs(verify['ra'] - angtracks[mjd][obj]['ra'][0]) < 60. / 3600]
                            for obj2 in verify['objid']:
                                if mjdalltracks[mjd][obj2]['ra'][0] in angtracks[mjd][obj]['ra'] or \
                                                mjdalltracks[mjd][obj2]['ra'][1] in angtracks[mjd][obj]['ra']:
                                    if set(mjdalltracks[mjd][obj2]['ra']) <= set(angtracks[mjd][obj]['ra']):
                                        ifsub.append(obj2)
                            if ifsub != []:
                                for delobj in ifsub:
                                    del mjdalltracks[mjd][delobj]
                                    keepbool = tracksarray['objid'] != delobj
                                    tracksarray = tracksarray[keepbool]
                            mjdalltracks[mjd][obj] = angtracks[mjd][obj]
                            tracksarray = np.append(tracksarray, np.array(
                                [(angtracks[mjd][obj]['ra'][0], angtracks[mjd][obj]['dec'][0], mjd, obj)],
                                dtype=tracksarray.dtype))
                        else:
                            mjdalltracks[mjd] = {obj: angtracks[mjd][obj]}
                            tracksarray = np.append(tracksarray, np.array(
                                [(angtracks[mjd][obj]['ra'][0], angtracks[mjd][obj]['dec'][0], mjd, obj)],
                                dtype=tracksarray.dtype))
            print '[%s] CleanTracks - Processing angle: %s. Total time:%s' % (
            time.strftime("%D %H:%M:%S"), n + 1, time.time() - t2)
        # print len(mjdalltracks[mjdalltracks.keys()[0]])
        self.mjdalltracks = mjdalltracks


def unwrap_self_f(arg, **kwarg):
    return SearchMovingObjects.oneangle(*arg, **kwarg)


def main():
    parser = argparse.ArgumentParser(
        description='Search a Non-Stationary DAOMOP HPX catalog for moving sources.')
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
    parser.add_argument("--rates-file", default="ratelist.17BQ02.txt",
                        help="File coontaining listing of rate and angle limits for search")
    parser.add_argument("qrunid", help="The CFHT QRUN to search.")

    cmd_line = " ".join(sys.argv)
    args = parser.parse_args()

    util.set_logger(args)
    logging.info("Started {}".format(cmd_line))

    storage.DBIMAGES = args.dbimages

    catalog_dirname = "{}/{}".format(args.catalogs, args.qrunid)
    search_catalog = storage.HPXCatalog(pixel=args.healpix, catalog_dir=catalog_dirname)
    trim_catalog = np.all((search_catalog.table['QRUNID'] == args.qrunid,
                           search_catalog.table['MATCHES'] == 0,
                           search_catalog.table['OVERLAPS'] > 1,
                           search_catalog.table['MAGERR_PSF'] < 0.25), axis=0)
    search_catalog.table = search_catalog.table[trim_catalog]
    search_catalog.table.sort("mid_mjdate")
    logging.info(len(search_catalog.table))

    #rates_listing = storage.Artifact(storage.Observation("rateslist"), version=args.qrunid, ext=".txt")
    #rates_listing.get()
    #rates_listing.table = ascii.read(rates_listing.filename)

    s = SearchMovingObjects(search_catalog)

    s.search()
    # s.CleanTracks_array()
    # json.dump(s.mjdalltracks, open(alltrackname, 'w'))

if __name__ == '__main__':
    main()
