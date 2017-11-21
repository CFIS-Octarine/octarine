import os
import sys
import time
import math
#from math import *
import copy
import json
import itertools
import subprocess
import warnings
warnings.filterwarnings("ignore")

#import pyfits
from astropy.io import fits as pyfits
import ephem
import numpy as np
import sqlite3 as sq3
from scipy.spatial import cKDTree as KDTree
#from contextlib import closing
from multiprocessing import Pool

__version__ = "2.4"

#####################################
######## Miscellaneous functions
#####################################

def arc(RA, DEC, pRA, pDEC):
	if RA == pRA and DEC == pDEC:
		return 0.0
	else:
		a = math.cos(math.radians(90.-DEC))*math.cos(math.radians(90.-pDEC))
		b = math.sin(math.radians(90.-DEC))*math.sin(math.radians(90.-pDEC)) * math.cos(math.radians(RA-pRA))
		if a + b > 1:
			return 0.0
		else:
			arc = math.degrees(math.acos(a+b))
			return arc

class OutMPCFormat:
	def __init__(self, mjd, ra, dec, mag, filterID, obscode):
		self.mjd = mjd
		self.ra = ra
		self.dec = dec
		self.mag = mag
		self.filterID = filterID
		self.obscode = obscode
	def mjd2date(self):
		date = ephem.date(self.mjd-15019.5)
		(year, mon, day) = date.triple()
		return "%s %02d %09.6f" %(year,mon,day)
	def make_radec(self, x, y):
		#make ra
		x = float(x)/15.
		ra = "%.2i %.2i %06.3f" %(int(x//1), x%1*60//1, x%1*60%1*60)
		#make dec
		y = float(y)
		if y > 0 :
			PN = "+"
		else:
			PN = "-"
		y = abs(y)
		dec = "%s%.2i %.2i %05.2f"  %(PN, int(y//1), y%1*60//1, y%1*60%1*60)
		return ra, dec
	def OutMPCline(self):
		date = self.mjd2date()
		s_ra, s_dec = self.make_radec(self.ra, self.dec)
		for n in range(1,5):
			if not self.filterID[-n].isdigit():
				band = self.filterID[-n]
				break
		out = "     NONE     C%s%s%s         %5.2f%s      %s\n" %(date, s_ra, s_dec, self.mag, band, self.obscode)
		return out

class OrbfitDiagnose:
	def __init__(self, mpclines, minau, maxau, crit_residual, printmpc):
		self.mpclines = mpclines
		self.minau = minau
		self.maxau = maxau
		self.crit_residual =crit_residual
		self.printmpc = printmpc
	def check_residule(self, stderr, crit_residual):
		for k in stderr:
			k = k.split()
			if k[0].isdigit():
				if abs(float(k[3])) > crit_residual or abs(float(k[5])) > crit_residual:	
					return False
		return True
	def Get_a_e_i(self, stdout):
		orbit = filter(lambda x: x.startswith("# a="), stdout)[0].strip()
		baryline = filter(lambda x: x.startswith("# Barycentric distance"), stdout)[0].strip()
		bary, baryerr =  float(baryline.split()[3].split('+-')[0]), float(baryline.split()[3].split('+-')[1])
		au	= float(orbit.split()[1].split("=")[1])
		e	= float(orbit.split("e=")[1].split(',')[0])
		inc	= float(orbit.split("i=")[1].split()[0])
		chi_line = filter(lambda x: x.startswith('# Chi-squared of fit:'), stdout)[0]
		chi, DOF = float(chi_line.split()[4]), float(chi_line.split()[6])
		rchi = chi/DOF
		return au, e, inc, chi, DOF, rchi, bary, baryerr
	def orbfit(self):
		fitout = subprocess.Popen("printf '%s' | fit_radec" %(self.mpclines), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = fitout.stdout.readlines(), fitout.stderr.readlines()
		fitgood = False
		if '# fit_radec\n' in stdout and "Best fit orbit gives:\n" in stderr:
			au, e, inc, chi, DOF, rchi, bary, baryerr = self.Get_a_e_i(stdout)
			if self.maxau > au > self.minau and self.check_residule(stderr, self.crit_residual):
				if self.printmpc == "yes":
					print "#%s %s %s %s %s %s %s %s" %(au, e, inc, chi, DOF, rchi, bary, baryerr)
					print self.mpclines.rstrip()
					for err in stderr:
						print err.strip()
					print ""
				fitgood = True
		else:
			print "===========fit_radec error !!!!!!==========="
			print "mpc lines:\n%sstdout:\n\t%s\nstderr:\n\t%s" %(self.mpclines, stdout, stderr)
			au, e, inc = 9999.9, 0.99999, 999.99
			pass
		return fitgood, au, e, inc



class linking:
	def __init__(self, mjdalltracks, healpix):
		self.mjdalltracks = mjdalltracks
		self.healpix = healpix
		self.processlist = list(itertools.combinations(self.mjdalltracks.keys(),2))
		self.run = "%s.r%3.1f.a%03i-%03i" %(self.healpix, crit_residual, int(maxau), int(minau))
	def FindCandidates(self):
		vslow = 147./slowau*24/3600.
		vslowabs = abs(vslow)
		t1= time.time()
		output = []
		processlist = list(itertools.combinations(self.mjdalltracks.keys(),2))
		for mm in processlist:
			mmm  = list(mm)
			mmm.sort()
			(mjd1, mjd2) = mmm[0], mmm[1]
			print "[%s]Process %i-%i [%s:%s]" %(time.strftime("%D %H:%M:%S"), int(mjd1), int(mjd2), len(self.mjdalltracks[mjd1].keys()), len(self.mjdalltracks[mjd2].keys()))
			oo=0
			for obj1name in self.mjdalltracks[mjd1].keys():
				oo +=1
				# print oo, len(self.mjdalltracks[mjd1].keys())
				oo2 = 0
				for obj2name in self.mjdalltracks[mjd2].keys():
					oo2+=1
					obj1 = self.mjdalltracks[mjd1][obj1name]
					obj2 = self.mjdalltracks[mjd2][obj2name]
					dt1 = (obj1['mjd'][-1] - obj1['mjd'][0])
					dt2 = (obj2['mjd'][-1] - obj2['mjd'][0])
					vobj1ra = ( obj1['ra'][-1] - obj1['ra'][0]) / dt1
					vobj2ra = ( obj2['ra'][-1] - obj2['ra'][0]) / dt2
					# only check the pairs with the same direction of vector
					if vobj1ra * vobj2ra > 0:
						dtobj12 = obj2['mjd'][0] - obj1['mjd'][0]
						predict_err = abs(dtobj12 / dt1) * astrometryerror
						vobj1dec = ( obj1['dec'][-1] - obj1['dec'][0]) / dt1
						vobj2dec = ( obj2['dec'][-1] - obj2['dec'][0]) / dt2
						x_pred = obj1['ra'][0] + vobj1ra * dtobj12
						y_pred = obj1['dec'][0] + vobj1dec * dtobj12
						vobj1 = (vobj1ra**2 + vobj1dec**2)**0.5
						vobj2 = (vobj2ra**2 + vobj2dec**2)**0.5
						# only try BK if second track within the astrometry error
						goodpair = False
						if ( abs(obj2['ra'][0] - x_pred)**2 + abs(obj2['dec'][0] - y_pred)**2)**0.5 < predict_err:
							if vobj1 < vslowabs and vobj2 < vslowabs:
								goodpair = True
							else:
								if abs((vobj2-vobj1)/min([vobj1, vobj2])) < 1 :
									goodpair = True
						if goodpair:
							out = ""
							obj1mpclist, obj2mpclist = [], []
							for n, x in enumerate(obj1['ra']):
								m = OutMPCFormat(obj1['mjd'][n], obj1['ra'][n], obj1['dec'][n], obj1['mag'][n], obj1['filterid'][n], obscode)
								ml = m.OutMPCline()
								out = out + ml
								obj1mpclist.append(ml)
							for n, x in enumerate(obj2['ra']):
								m = OutMPCFormat(obj2['mjd'][n], obj2['ra'][n], obj2['dec'][n], obj2['mag'][n], obj2['filterid'][n], obscode)
								ml = m.OutMPCline()
								out = out + ml
								obj2mpclist.append(ml)
							c_orbfit = OrbfitDiagnose(out, minau, maxau, crit_residual, printmpc)
							fitgood, a, ecc, inc = c_orbfit.orbfit()
							if fitgood:
								print oo, oo2, obj1name, obj2name
								saveinfo = []
								for n, x in enumerate(obj1['ra']):
									saveinfo.append([obj1['ra'][n], obj1['dec'][n], obj1['mjd'][n], obj1['mag'][n], obj1['magerr'][n], obj1['filterid'][n], obscode, obj1['exptime'][n], obj1['fitsname'][n], a, ecc, inc, obj1mpclist[n]])
								for n, x in enumerate(obj2['ra']):
									saveinfo.append([obj2['ra'][n], obj2['dec'][n], obj2['mjd'][n], obj2['mag'][n], obj2['magerr'][n], obj2['filterid'][n], obscode, obj2['exptime'][n], obj2['fitsname'][n], a, ecc, inc, obj2mpclist[n]])
								output.append(saveinfo)
							else:
								output.append([])
		self.output = output
	def savedb(self):		
		candnum = 1
		db = MovingObjectDB(self.run)
		for i in self.output:
			if i !=[]:
				cand = '%s_cand%05d' %(self.healpix, candnum)
				print "%s, a=%s, e=%s, i=%s, mag=%s" %(cand, i[0][9], i[0][10], i[0][11], i[0][3])
				for j in i:
					[x1, y1, mjd1, mag1, magerr1, band1, obscode, exp1, chip1, a, ecc, inc, out1] = j
					db.insertdb(cand, x1, y1, mjd1, mag1, magerr1, band1, obscode, exp1, chip1, a, ecc, inc, out1)
				candnum +=1
# 				[x1, y1, mjd1, mag1, magerr1, band1, obscode, exp1, chip1, a, ecc, inc, out1] = i[0]
# 				[x2, y2, mjd2, mag2, magerr2, band2, obscode, exp2, chip2, a, ecc, inc, out2] = i[1]
# 				[x3, y3, mjd3, mag3, magerr3, band3, obscode, exp3, chip3, a, ecc, inc, out3] = i[2]
# 				[x4, y4, mjd4, mag4, magerr4, band4, obscode, exp4, chip4, a, ecc, inc, out4] = i[3]
# 				db.insertdb(cand, x1, y1, mjd1, mag1, magerr1, band1, obscode, exp1, chip1, a, ecc, inc, out1)
# 				db.insertdb(cand, x2, y2, mjd2, mag2, magerr2, band2, obscode, exp2, chip2, a, ecc, inc, out2)
# 				db.insertdb(cand, x3, y3, mjd3, mag3, magerr3, band3, obscode, exp3, chip3, a, ecc, inc, out3)
# 				db.insertdb(cand, x4, y4, mjd4, mag4, magerr4, band4, obscode, exp4, chip4, a, ecc, inc, out4)		
# 				candnum +=1
		db.close()
		#print "%s, %s, %s" %(self.ori, self.cdeg, self.of)


#####################################
### The class for generating synthetic SS objects
#####################################

class make_planet:
	def __init__(self, ra, dec, mjd, a0, e0, i0, H):
		self.a0 = a0
		self.e0 = e0
		self.i0 = i0
		self.H = H
		self.mjd = mjd
		self.ra = ra
		self.dec = dec
		np = ephem.Equatorial(math.radians(ra), math.radians(dec), epoch='2000')
		e = ephem.Ecliptic(np)
		self.e_ra = math.degrees(e.lon)
		self.e_dec = math.degrees(e.lat)
	def planet(self):
		name = 'object' 
		a = self.a0
		e = self.e0
		i = self.i0
		found = False
		for s in np.arange(-180, 180, 1):
			if self.e_dec > 0:
				O = self.ra - self.e_dec / (math.tan(math.radians(i))) + s 
			else:
				O = self.ra + self.e_dec / (math.tan(math.radians(i))) + s 
			for p in np.arange(-180, 180, 1):
				o = self.e_dec/math.sin(math.radians(i)) +p
				M = 0.01
				H = self.H
				n = 0
				# JD = self.jd
				d = ephem.Date(self.mjd-15019.5)
				E = '%i/%i/%i' %(int(d.triple()[1]), d.triple()[2], d.triple()[0])
				D = '2000.0'
				G = '0.15'
				output = "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" %(name, 'e', i, O, o, a, n, e, M, E, D, H, G)
				orbit = ephem.readdb("'%s'" %(output))
				orbit.compute(epoch='2000')
				orbit.compute(self.mjd - 15019.5)
				mag = float(orbit.mag)
				ora = math.degrees(orbit.a_ra)
				odec = math.degrees(orbit.a_dec)
#				print ora, odec, O, o
				if abs(ora-self.ra) < 0.5 and abs(odec - self.dec) < 0.5:
					break
			if abs(ora-self.ra) < 0.5 and abs(odec - self.dec) < 0.5:
				found = True
				break
		return output, found


#####################################
### Build db
#####################################
class MovingObjectDB:
	def __init__(self,run):
		self.dbdir = os.getcwd()
		self.namedb = 'HSC_SS.%s.db' %(run)
		self.nametb = 'cands'
		self.dbpath = "%s/%s" %(self.dbdir, self.namedb)
		self.creatDB()
		self.conn = sq3.connect(self.dbpath)
		self.cursor = self.conn.cursor()
	def creatDB(self):
		if os.path.exists(self.dbpath):
			os.system('rm %s' %(self.dbpath))
		conn = sq3.connect(self.dbpath)
		cursor = conn.cursor()
		cursor.execute ("""
CREATE TABLE %s(cand text, 
  RA  real,
  DEC  real,
  MJD  real,
  mag  real,
  magerr  real,
  filterid  text,
  obscode  text,
  exptime real,
  frameid  text,
  a  real,
  ecc  real,
  inc  real,
  mpcline  text)""" %(self.nametb))
		cursor.close()
	def insertdb(self, cands, ra, dec, mjd, mag, magerr, filterid, obscode, exptime, frameid, a, ecc, inc, mpcline):
		self.cursor.execute("insert into %s values ('%s', %f, %f, %f, %f, %f, '%s', '%s', %f, '%s', %f, %f, %f, '%s')" %(self.nametb, cands, ra, dec, mjd, mag, magerr, filterid, obscode, exptime, frameid, a, ecc, inc, mpcline))
		self.conn.commit()
	def close(self):
		self.cursor.close()


#####################################
### Searching moving objects
#####################################

#HPX_02937_RA_160.3_DEC_+31.4_cat.fits.ns

class SearchMovingObjects:
	def __init__(self, workdir, healpix):
		self.workdir = workdir
		self.healpix = healpix
		os.chdir(self.workdir)
		self.cra = float(self.healpix.split('_')[0])
		self.cdec = float(self.healpix.split('_')[1])
	def CombineNS(self):
		# workdir = '/array/users/ytchen/processed/2016A/201606/rerun/ssp-201606-all-lownoise/healpix/345.9_+001.2_04711'
		print '[%s] Combining reasonable ns files... (%s)' %(time.strftime("%D %H:%M:%S"), self.workdir)
		# os.chdir(workdir)
		# os.system('rm %s' %(alldetname))
		ns2list = filter(lambda x: x.startswith('det') and x.endswith('fits') and 'mlns' in x, os.listdir(self.workdir))
		ns2list.sort()
		#if os.path.exists(alldetname):
		#	os.remove(alldetname)
		#	print '[%s] %s exists, deleting that...' %(time.strftime("%D %H:%M:%S"), alldetname)
		#os.system('touch %s' %(alldetname))
		jdlist = []
		for ns2name in ns2list:
			ptname = '%spointing' %(ns2name.rstrip('.fits').rstrip('mlns'))
			pt = open(ptname).readline()
			# 2457554.12075 200.0 HSC-Y 344.417683333 1.31750277778 2016/6/14 14:53:53 260.986707443
			jd = pt.split()[0]
			exptime = float(pt.split()[1])
			filterid = pt.split()[2]
			jdlist.append(float(jd))
			# Only choose gri bands and exptime > 100 seconds
			if filterid in ProcessFilterid and exptime > ProcessMinExptime and ifmakenewstackfits:
				ft = pyfits.open(ns2name)[1].data
				if 'outarray' in locals():
					outarray = np.concatenate(( outarray, ft ),axis=0)
				else:
					outarray = ft
			else:
	#			print 'passing', pt
				pass
		jdlist.sort()
		midjd = jdlist[int(len(jdlist)/2.)]
		self.alldets = outarray
		return midjd
	def readalldet(self):
		#sns = np.loadtxt(alldetname, dtype={'names':('ra', 'dec', 'jd', 'mag', 'magerr','A','B', 'theta', 'id', 'parent', 'fitsname', 'filterid', 'exptime'), 'formats':('f8','f8','f8', 'f8','f8','f8','f8','f8', '>i8', '>i8', 'S150', 'S10', 'f8')})
		# This condition pick mag > 26.0 and include 'nan' mag
		# bigmag = sns['mag'] >26.0
		# sns = sns[~bigmag]
		# This condition Only pick the reasonable mag < 26.0 (no 'nan' mag)
		m = -2.5*np.log10(self.alldets['flux_psf']/self.alldets['fz0']) 
		self.sns = self.alldets[m < 26.0]
	def readns(self, nsfile):
		ft = pyfits.open(nsfile)[1].data
		mjdlist = list(set([int(i) for i in ft['mid_mjdate']]))
		mjdlist.sort()
		midmjd = mjdlist[int(len(mjdlist)/2.)]
		outarray = ft[(ft['FLUX_RADIUS'] > 2)&(ft['MAG_ISO'] <24.5)&(ft['X_IMAGE'] < 2085) &(ft['X_IMAGE'] > 42)]
		self.sns = outarray
		global sns
		sns = ft[(ft['FLUX_RADIUS'] > 2)&(ft['MAG_ISO'] <24.5)&    (ft['X_IMAGE'] < 2085) &(ft['X_IMAGE'] > 42)]
		time.sleep(1)
		print 'self.sns, ft', len(self.sns), len(ft), len(outarray)
		return midmjd
	def FindMovingAngle(self, mjd, ra, dec):
		print '[%s] Getting the reasonable moving speed on (%s, %s) at %s (%s)' %(time.strftime("%D %H:%M:%S"), ra, dec, ephem.Date(mjd-15019.5), mjd)
		anglelist = []
		fp = ephem.Equatorial(math.radians(ra), math.radians(dec), epoch='2000')
		ep = ephem.Ecliptic(fp)
		e_ra = math.degrees(ep.lon)
		e_dec = math.degrees(ep.lat)
		for au in np.arange(au_s, au_e, 1.0):
			#for i in np.arange(e_dec, i_e, 2.0):
			for i in np.arange(e_dec, e_dec+10, 2.0):
				a = make_planet(ra, dec, mjd, au, e_s, i, H)
				l, found = a.planet()
				if found:
					orbit = ephem.readdb("'%s'" %(l))
					orbit.compute(epoch='2000')
					orbit.compute(mjd-15019.5)
					mag = float(orbit.mag)
					ra0 = math.degrees(orbit.a_ra)
					dec0 = math.degrees(orbit.a_dec)
					orbit.compute(mjd-15019.5+0.1)
					mag = float(orbit.mag)
					ra1 = math.degrees(orbit.a_ra)
					dec1 = math.degrees(orbit.a_dec)
					ang = math.degrees(math.atan((dec1-dec0)/(ra1-ra0)))
					anglelist.append(ang)
					#print au, i, ang
		hist, bin_edges = np.histogram(anglelist, bins=np.arange(-60, 60, 2))
		aveang = int((bin_edges[np.where(hist == hist.max())[0][0]] + bin_edges[np.where(hist == hist.max())[0][0]+1] ) /2.)
		print '[%s] Getting the reasonable moving speed on (%s %s)... Done. Mean angle = %s' %(time.strftime("%D %H:%M:%S"), ra, dec, aveang)
	#	return anglelist
		return aveang
	def searchline(self, inputra, inputdec, ang, sr):
		t2 = time.time()
		slope = np.tan(np.radians(ang))
		b = inputdec - slope * inputra
		snsra, snsdec = self.sns['X_WORLD'], self.sns['Y_WORLD']
		inline = abs(slope*snsra - snsdec + b)/((slope**2+1)**0.5) < sr/3600
		b = sns[ inline]
		if inline.sum()==0:
			return 0
		else:
			#print 'sns, self.sns, inline:', len(sns), len(self.sns), len(inline)
			#match = self.sns[inline]
			match = sns[inline]
			return match
	def oneangle(self, ang):
		t2 = time.time()
		result2 = []
		print '[%s] Processing angle: %s' %(time.strftime("%D %H:%M:%S"), ang)
		for sl, ldec in enumerate(np.arange(self.cdec - fr, self.cdec + fr, step/3600)):
			t3 = time.time()
			inputra, inputdec = self.cra, ldec
			mat = self.searchline(inputra, inputdec, ang, sr)
			if mat !=0:
				#print mat
				pairs = self.findpair(sl, mat, ang, fastau, slowau)
				result2.append(pairs)
			#print time.time()-t3
		print '[%s] Processing angle: %s. Total pair: %s. Total time: %s' %(time.strftime("%D %H:%M:%S"), ang, len(result2), time.time()-t2)
		return result2
	def findpair(self,sl, mat, ang, fastau, slowau):
		vmax = 147./fastau*24/3600
		vmaxabs = abs(vmax)
		vslow = 147./slowau*24/3600.
		vslowabs = abs(vslow)
		line = mat
		#mjdlist = list(set([int(i - 2400000.5) for i in line['jd']]))
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
			for n, i in enumerate(sorted(mdata, key = lambda mdata: np.degrees(mdata['X_WORLD']))):
				#tra, tdec, tjd, tmag, tmagerr = i['ra'], i['dec'], i['jd'], i['mag'], i['magerr'] 
				#tA, tB, ttheta, tpath = i['A'], i['B'], i['theta'], i['fitsname']
				#tf, tet = i['filterid'], i['exptime']
				#tra, tdec, tmjd = i['X_WORLD'], i['Y_WORLD'], i['mid_mjdate']
				tra, tdec, tmjd = float(i['X_WORLD']), float(i['Y_WORLD']), float(i['mid_mjdate'])
				#tmag = -2.5*np.log10(i['flux_psf']/i['fz0']) 
				#tmagerr = -2.5*np.log10(1-i['flux_psf_err']/i['flux_psf'])
				#tmag = i[Processmag]
				#tmagerr = i[Processmagerr]
				tmag = float(i[Processmag])
				tmagerr = float(i[Processmagerr])
				#shape = i['shape_sdss'] # unit: pixel^2
				#xx, yy, xy = shape[0], shape[1], shape[2]
				#tA = (((xx+yy)/2) + (  ((xx-yy)/2)**2 + xy**2)**0.5)**0.5
				#tB = (((xx+yy)/2) - (  ((xx-yy)/2)**2 + xy**2)**0.5)**0.5
				#ttheta = np.degrees(1./2*np.arctan2((2*xy),(xx-yy)))
				#tA = i['A_IMAGE']
				#tB = i['B_IMAGE']
				#ttheta = i['THETA_IMAGE']
				tA = float(i['A_IMAGE'])
				tB = float(i['B_IMAGE'])
				ttheta = float(i['THETA_IMAGE'])
				tpath, tf, tet = i['dataset_name'], 'r', 30.0
				radecstr = '%13.9f_%13.11f_%03i_%04i_%04i' %(tra, tdec, ang, sl, n)
				linking[radecstr] = {'ra':[tra], 'dec':[tdec], 'mjd':[tmjd], 'mag':[tmag], 'magerr':[tmagerr], 'A':[tA], 'B':[tB], 'theta':[ttheta], 'fitsname':[tpath], 'filterid':[tf], 'exptime':[tet]}
				for obj in linking.keys():
					tk = linking[obj]
					if tmjd != tk['mjd'][0] and tmjd != tk['mjd'][-1]:
						a_dt = abs(tmjd - tk['mjd'][-1])
						d2 = ((tdec - tk['dec'][0])**2 + (tra - tk['ra'][0])**2)**0.5
						a_v = d2 / a_dt
						if a_dt > min_dt and a_v < max_v and a_v > min_v:
							motion = ((tdec - tk['dec'][-1])**2 + (tra - tk['ra'][-1])**2)**0.5
							if motion > min_motion/3600.:
								# going to add detection after detection, no matter prograde or retrograde
								v2 = d2/(tmjd - tk['mjd'][0])
								v2abs = abs(v2)
								newobj = '%s_%04i' %(obj, n)
								if abs(v2) < abs(vmax) and abs(tk['mag'][0] - tmag) < max_dm:
									if tk['ra'].__len__() ==1:
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
										#moving = False
										#if tk['jd'][-1] - tk['jd'][0] > 0 and tjd > tk['jd'][-1]:
										#	moving = True
										#if tk['jd'][-1] - tk['jd'][0] < 0 and tjd < tk['jd'][-1]:
										#	moving = True
										#if moving:
										# if the direction of motions are the same in RA
										if (tk['mjd'][-1] - tk['mjd'][0]) * (tmjd - tk['mjd'][-1]) > 0:
											goodvelocity = False
											# The condition for the larger error ratio of slow moving objects in short duration
											d1 = ((tk['dec'][-1] - tk['dec'][0])**2 + (tk['ra'][-1] - tk['ra'][0])**2 )**0.5
											v1 = d1/(tk['mjd'][-1] - tk['mjd'][0])
											v1abs = abs(v1)
											# condition for slow mover (intra-night)
											if v1abs < vslowabs and v2abs < vslowabs:
												#if abs( (v2-v1)/min([v1, v2]) ) < 1 and abs(v1) < vmaxabs:
												goodvelocity = True
											# condition for fast mover (intra-night)
											else:
												if abs( (v2-v1)/min([v1abs, v2abs]) ) < 0.2 and v1abs < vmaxabs:
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
			#linkingout = copy.deepcopy(linking)
			if deduplicate:
				dlist = []
				for obj1 in linking.keys():
					if linking[obj1]['ra'].__len__() < N_dets or max(linking[obj1]['mag']) - min(linking[obj1]['mag']) > max_dm:
						#del linkingout[obj1]
						dlist.append(obj1)
				for obj1 in dlist:
					del linking[obj1]
				#linking = copy.deepcopy(linkingout)
				dlist = []
				for obj1 in linking.keys():
					for obj2 in linking.keys():
						if obj1 != obj2:
							if set(linking[obj1]['ra']) < set(linking[obj2]['ra']):
								#del linkingout[obj1]
								dlist.append(obj1)
								break
				for obj1 in dlist:
					del linking[obj1]				
			# pairs.append(linkingout)
			#pairs[m] = linkingout
			pairs[m] = linking
		return pairs
	def Searching(self, multi, nsfile):
		#self.midjd = self.CombineNS(True)
		#self.readalldet()
		self.midmjd = self.readns(nsfile)
		#self.aveang = self.FindMovingAngle(self.midmjd, self.cra, self.cdec)
#		self.aveang = 31
		# ratelist = open('/sciproc/disk2/cfis/mis/ratelist.txt').readlines()
		# ratelist = open('/sciproc/disk2/cfis/mis/ratelist.17BQ02.txt').readlines()
		ratelist = open(sys.argv[3]).readlines()
		d = {}
		for r in ratelist:
			rns = r.split()[0]
			rate =int(float(r.split()[2]))
			d[rns] = rate
		self.aveang = d[nsfile.rstrip('ns').rstrip('.')]
		print '[%s] Getting the reasonable moving speed on (%s)... Done. Mean angle = %s' %(time.strftime("%D %H:%M:%S"), nsfile, self.aveang)
		allanglist = np.arange(self.aveang-openangle, self.aveang+openangle, 1.0)
		pool = Pool(processes = N_threading)
		#self.subresults = pool.map(oneangle, allanglist)
		print '[%s] Searching reasonable tracks ... ' %(time.strftime("%D %H:%M:%S"))
		if multi:
			self.results = pool.map(unwrap_self_f, zip([self]*len(allanglist), allanglist), chunksize=1)
		else:
			self.results = []
			for ang in allanglist:
				r = self.oneangle(ang)
				self.results.append(r)	
	def CleanTracks(self):
		mjdalltracks = {}
		for n, angtrackslist in enumerate(s.results):
			t2 = time.time()
			for angtracks in angtrackslist:
				for mjd in angtracks.keys():
					for obj in angtracks[mjd]:
						# Checking if mjd is in mjdalltracks dictionary
						if mjd in mjdalltracks.keys():
							# Checking the deplicates
							ifsub = []
							for obj2 in mjdalltracks[mjd]:
								if mjdalltracks[mjd][obj2]['ra'][0] in angtracks[mjd][obj]['ra'] or mjdalltracks[mjd][obj2]['ra'][1] in angtracks[mjd][obj]['ra']:
									if set(mjdalltracks[mjd][obj2]['ra']) <= set(angtracks[mjd][obj]['ra']):
										ifsub.append(obj2)
							if ifsub !=[]:
								for obj2 in ifsub:
									del mjdalltracks[mjd][obj2]
								mjdalltracks[mjd][obj] = copy.deepcopy(angtracks[mjd][obj])
							else:
								mjdalltracks[mjd][obj] = copy.deepcopy(angtracks[mjd][obj])
						else:
							mjdalltracks[mjd] = {obj:copy.deepcopy(angtracks[mjd][obj])}
			print '[%s] CleanTracks - Processing angle: %s. Total time:%s' %(time.strftime("%D %H:%M:%S"), n+1, time.time()-t2)
		self.mjdalltracks = mjdalltracks
	def CleanTracks_array(self):
		mjdalltracks = {}
		tracksarray = np.array([ ], dtype = [('ra', 'f8'), ('dec', 'f8'), ('mjd', 'i4'), ('objid', 'S150')])
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
							verify = verify[abs(verify['ra'] - angtracks[mjd][obj]['ra'][0]) < 60./3600]
							for obj2 in verify['objid']:
								if mjdalltracks[mjd][obj2]['ra'][0] in angtracks[mjd][obj]['ra'] or mjdalltracks[mjd][obj2]['ra'][1] in angtracks[mjd][obj]['ra']:
									if set(mjdalltracks[mjd][obj2]['ra']) <= set(angtracks[mjd][obj]['ra']):
										ifsub.append(obj2)
							if ifsub !=[]:
								for delobj in ifsub:
									del mjdalltracks[mjd][delobj]
									keepbool = tracksarray['objid'] != delobj
									tracksarray = tracksarray[keepbool]
							mjdalltracks[mjd][obj] = angtracks[mjd][obj]
							tracksarray = np.append(tracksarray, np.array([(angtracks[mjd][obj]['ra'][0], angtracks[mjd][obj]['dec'][0], mjd, obj)], dtype=tracksarray.dtype))
						else:
							mjdalltracks[mjd] = {obj:angtracks[mjd][obj]}
							tracksarray = np.append(tracksarray, np.array([(angtracks[mjd][obj]['ra'][0], angtracks[mjd][obj]['dec'][0], mjd, obj)], dtype=tracksarray.dtype))
			print '[%s] CleanTracks - Processing angle: %s. Total time:%s' %(time.strftime("%D %H:%M:%S"), n+1, time.time()-t2)
			#print len(mjdalltracks[mjdalltracks.keys()[0]])
		self.mjdalltracks = mjdalltracks


def unwrap_self_f(arg, **kwarg):
    return SearchMovingObjects.oneangle(*arg, **kwarg)



#alldetname = 'alldet.ns2'

# parameters for generate synthetic SS objects
au_s = 3.0 # starting au
au_e = 40.0 # ending au
# i_e = 30.0 # ending inclination  # not use anymore
e_s = 0.11 # starting eccentricity
H = 16.0 # assuming SS absolute magintude

# search conditions for stacked ns catalog
openangle = 30.
fr = 2.5 # the field/healpix radius, in degree [2.5 for moving rate up to 75 degree, 2 for 60 degree, 1.3 for < 45]
sr = 5.0 # the search radius for searchline, in arcsecond
step = 2.0 # the search step for searchline, in arcsecond
fastau = 3.0 # the fastet moving speed for the object at * au, the tracks faster than this value will be ignored. 
slowau = 30.0 # the slow moving speed for the object at * au, the tracks slower than this value will be *USING LOOSER CONDITION OF SPEED*	to avoid over-kill.
N_dets = 3 # only consider the tracks with N_dets at least
max_dm = 1.0 # only consider the tracks with dm < max_dm
min_motion = 0.2 # in arcsecond, the requirement of minmum motion between two exposures 
min_dt  = 20./60/24 # in day, the requirement of minmum time separation.
max_v = 15.0  *24/3600 # degree per day
min_v = 1.0  *24/3600 # degree per day
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
astrometryerror = 1./3600 #arcsec

#ProcessFilterid = ['HSC-g', 'HSC-r', 'HSC-i', 'HSC-r2', 'HSC-i2']
Processmag = 'MAG_ISO'
Processmag = 'MAG_PSF'
Processmagerr = 'MAGERR_ISO'
Processmagerr = 'MAGERR_PSF'
ProcessMinExptime = 10.
stackfitsname = 'allmlns.fits'
allpairname = 'allpairs.json'
alltrackname = 'mjdalltracks.json'
scn = 'stationary_catalog_ndet_2_r0.5_dt20m'

if __name__ == '__main__':
	workdir = sys.argv[1]
	nsfile = sys.argv[2]
	alltrackname = '%s%s' %(nsfile.rstrip('ns'), alltrackname)
	#nsfile = 'HPX_02937_RA_160.3_DEC_+31.4'
	healpix = '%s_%s_%s' %(nsfile.split('_')[3], nsfile.split('_')[5], nsfile.split('_')[1] )
	os.chdir(workdir)
	s = SearchMovingObjects(workdir, healpix)
	if os.path.exists(allpairname) and os.path.exists(alltrackname):
		print '[%s] %s has been done. Loading the file' %(time.strftime("%D %H:%M:%S"), alltrackname)
		#s.results = json.load(open(allpairname, 'r'))
		s.mjdalltracks = json.load(open(alltrackname, 'r'))
	else:
		s.Searching(multi, nsfile)
		#json.dump(s.results, open(allpairname, 'w'))
		s.CleanTracks_array()
		json.dump(s.mjdalltracks, open(alltrackname, 'w'))
#	l = linking(s.mjdalltracks, healpix)
#	l.FindCandidates()
#	l.savedb()
	







