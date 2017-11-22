import os
import sys
import time
import copy
import json
import warnings
from astropy.io import fits as pyfits
import numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")

__version__ = "2.4"

# parameters for generate synthetic SS objects
au_s = 3.0  # starting au
au_e = 40.0  # ending au
# i_e = 30.0 # ending inclination  # not use anymore
e_s = 0.11  # starting eccentricity
H = 16.0  # assuming SS absolute magintude

# search conditions for stacked ns catalog
openangle = 30.
fr = 2.5  # the field/healpix radius, in degree [2.5 for moving rate up to 75 degree, 2 for 60 degree, 1.3 for < 45]
sr = 5.0  # the search radius for searchline, in arcsecond
step = 2.0  # the search step for searchline, in arcsecond
fastau = 3.0  # the fastet moving speed for the object at * au, the tracks faster than this value will be ignored.
slowau = 30.0  # the slow moving speed for the object at * au, slower tracks *LOOSen CONDITION OF SPEED *
N_dets = 3  # only consider the tracks with N_dets at least
max_dm = 1.0  # only consider the tracks with dm < max_dm
min_motion = 0.2  # in arcsecond, the requirement of minmum motion between two exposures
min_dt = 20. / 60 / 24  # in day, the requirement of minmum time separation.
max_v = 15.0 * 24 / 3600  # degree per day
min_v = 1.0 * 24 / 3600  # degree per day
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

# HPX_02937_RA_160.3_DEC_+31.4_cat.fits.ns

class SearchMovingObjects:
    def __init__(self, workdir, healpix):
        self.workdir = workdir
        self.healpix = healpix
        os.chdir(self.workdir)
        self.cra = float(self.healpix.split('_')[0])
        self.cdec = float(self.healpix.split('_')[1])

    def readns(self, nsfile):
        ft = pyfits.open(nsfile)[1].data
        mjdlist = list(set([int(i) for i in ft['mid_mjdate']]))
        mjdlist.sort()
        midmjd = mjdlist[int(len(mjdlist) / 2.)]
        outarray = ft[(ft['FLUX_RADIUS'] > 2) & (ft['MAG_ISO'] < 24.5) & (ft['X_IMAGE'] < 2085) & (ft['X_IMAGE'] > 42)]
        self.sns = outarray
        global sns
        sns = ft[(ft['FLUX_RADIUS'] > 2) & (ft['MAG_ISO'] < 24.5) & (ft['X_IMAGE'] < 2085) & (ft['X_IMAGE'] > 42)]
        time.sleep(1)
        print 'self.sns, ft', len(self.sns), len(ft), len(outarray)
        return midmjd

    def searchline(self, inputra, inputdec, ang, sr):

        slope = np.tan(np.radians(ang))
        b = inputdec - slope * inputra

        snsra, snsdec = self.sns['X_WORLD'], self.sns['Y_WORLD']
        inline = abs(slope * snsra - snsdec + b) / ((slope ** 2 + 1) ** 0.5) < sr / 3600.0

        if sum(inline) == 0:
            return 0
        else:
            match = sns[inline]
            return match

    def oneangle(self, ang):
        t2 = time.time()
        result2 = []
        print '[%s] Processing angle: %s' % (time.strftime("%D %H:%M:%S"), ang)
        for sl, ldec in enumerate(np.arange(self.cdec - fr, self.cdec + fr, step / 3600)):
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
                                            # The condition for the larger error ratio of slow moving objects in short duration
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

    def Searching(self, multi, nsfile):
        # self.midjd = self.CombineNS(True)
        # self.readalldet()
        self.midmjd = self.readns(nsfile)
        # self.aveang = self.FindMovingAngle(self.midmjd, self.cra, self.cdec)
        #		self.aveang = 31
        # ratelist = open('/sciproc/disk2/cfis/mis/ratelist.txt').readlines()
        # ratelist = open('/sciproc/disk2/cfis/mis/ratelist.17BQ02.txt').readlines()
        ratelist = open(sys.argv[3]).readlines()
        d = {}
        for r in ratelist:
            rns = r.split()[0]
            rate = int(float(r.split()[2]))
            d[rns] = rate
        self.aveang = d[nsfile.rstrip('ns').rstrip('.')]
        print '[%s] Getting the reasonable moving speed on (%s)... Done. Mean angle = %s' % (
        time.strftime("%D %H:%M:%S"), nsfile, self.aveang)
        allanglist = np.arange(self.aveang - openangle, self.aveang + openangle, 1.0)
        pool = Pool(processes=N_threading)
        # self.subresults = pool.map(oneangle, allanglist)
        print '[%s] Searching reasonable tracks ... ' % (time.strftime("%D %H:%M:%S"))
        if multi:
            self.results = pool.map(unwrap_self_f, zip([self] * len(allanglist), allanglist), chunksize=1)
        else:
            self.results = []
            for ang in allanglist:
                r = self.oneangle(ang)
                self.results.append(r)

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
    workdir = sys.argv[1]
    nsfile = sys.argv[2]
    global alltrackname
    alltrackname = '%s%s' % (nsfile.rstrip('ns'), alltrackname)
    # nsfile = 'HPX_02937_RA_160.3_DEC_+31.4'
    healpix = '%s_%s_%s' % (nsfile.split('_')[3], nsfile.split('_')[5], nsfile.split('_')[1])
    os.chdir(workdir)
    global s
    s = SearchMovingObjects(workdir, healpix)
    if os.path.exists(allpairname) and os.path.exists(alltrackname):
        print '[%s] %s has been done. Loading the file' % (time.strftime("%D %H:%M:%S"), alltrackname)
        # s.results = json.load(open(allpairname, 'r'))
        s.mjdalltracks = json.load(open(alltrackname, 'r'))
    else:
        s.Searching(multi, nsfile)
        # json.dump(s.results, open(allpairname, 'w'))
        s.CleanTracks_array()
        json.dump(s.mjdalltracks, open(alltrackname, 'w'))

if __name__ == '__main__':
    main()
