from astropy.io import fits
from matplotlib import pyplot as plt
from numpy import *
from sklearn.ensemble import RandomForestClassifier as RF
from sklearn.ensemble import IsolationForest as IF
from sklearn.metrics import confusion_matrix, roc_curve
from sklearn.model_selection import train_test_split
from scipy.spatial import cKDTree as KDTree
import pandas as pd
import seaborn as sns
import os, sys

class fits_cat:
    def __init__(self, fitsfile):
        self.fits_name = fitsfile
        self.fits_file = fits.open(self.fits_name)
        self.header = self.fits_file[1].header
        self.catalog = self.fits_file[1].data
        self.stationary0 = self.catalog[self.catalog['MATCHES'] > 0]
        self.stationary1 = self.catalog[self.catalog['MATCHES'] > 1]
        self.stationary2 = self.catalog[self.catalog['MATCHES'] > 2]
        self.ns0 = self.catalog[logical_and(logical_and(self.catalog['MATCHES'] == 0, self.catalog['OVERLAPS'] >0), self.catalog['MAG_PSF']<27)]
        self.ns1 = self.catalog[logical_and(self.catalog['MATCHES'] == 0, self.catalog['OVERLAPS'] >1)]
        self.ns2 = self.catalog[logical_and(self.catalog['MATCHES'] == 0, self.catalog['OVERLAPS'] >2)]
	self.real_feature = self.feature(self.stationary2)
        self.ns_feature = self.feature(self.ns0)
    def plot(self):
        #plt.xlim(165.5, 168.0)
        #plt.ylim(31.5, 32.0)
        #plt.xlim(165, 168.0)
        #plt.plot(self.stationary0['X_WORLD'], self.stationary0['Y_WORLD'], '.', label='MATCH >0')
        #plt.plot(self.stationary1['X_WORLD'], self.stationary1['Y_WORLD'], '.', label='MATCH >1')
        #plt.plot(self.ns['X_WORLD'], self.ns['Y_WORLD'], '.', label='NS, OVERLAPS >1')
        #plt.plot(self.ns1['X_WORLD'], self.ns1['Y_WORLD'], '.', label='NS, OVERLAPS >2')
        plt.plot(self.catalog['X_WORLD'], self.catalog['Y_WORLD'], marker = '.', label='ALL')
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
        
class IF_real_bogus:
    def __init__(self, feature):
        self.feature = feature
        #self.Class = []
        #[self.Class.append('real') for i in range(len(self.feature))]
        #self.trainning, self.test, self.Class_train, self.Class_test = train_test_split(self.real_feature,\
        #                        self.real_Class, test_size=0.4, random_state=0)
        #print self.trainning     
    def train(self):
        self.IFmod = IF(n_estimators = 160)
        self.IFmod.fit(self.feature)
        
    def validation(self):
        #score = self.IFmod.decision_function(self.feature)
        result = self.IFmod.predict(self.feature)
        self.normal = (result == 1)
        self.abnormal = (result == -1)
        print "total:{0}, normal:{1}, abnormal:{2}, normal/total:{3}".format(len(self.feature), 
            	self.normal.sum(), self.abnormal.sum(), self.normal.sum()/float(len(self.feature)))
                    
def build_model(fits_tab):
	sources = fits_cat(fits_tab)
	#print len(sources.stationary2), len(sources.ns1)
	real_bogus = IF_real_bogus(sources.real_feature)
	real_bogus.train()
	real_bogus.validation()
	return real_bogus.IFmod

def remove_ps1_source(ps1_cat, ns_cat):
	ps1 = fits.open(ps1_cat)[1].data
	ps1_ra = ps1['f0']
	ps1_dec = ps1['f1']
	ra = ns_cat['X_WORLD']
	dec = ns_cat['Y_WORLD']
	ns_tree = KDTree(zip(ra, dec))
	ps1_tree = KDTree(zip(ps1_ra, ps1_dec))
	match = ns_tree.query_ball_tree(ps1_tree, r=.5/3600.)
	#print len(match), len(ra), len(ps1)
	ps1_match = []
	for i in match:
		if i != []:
			ps1_match.append(True)
		else:
			ps1_match.append(False)
	return array(ps1_match)
	

def ns_cat(fits_tab):
	try:
		sources = fits_cat(fits_tab)
		healpix = fits_tab.split('/')[-1].rstrip('_cat.fits') #HPX_03205_RA_194.1_DEC_+28.6_cat.fits
		print fits_tab
		real_bogus = IF_real_bogus(sources.ns_feature)
		real_bogus.IFmod = RB_model
		real_bogus.validation()
		ns_cat = sources.ns0[real_bogus.normal]
		ps1_cat = '{0}{1}_STK.fits'.format(ps1_dir, healpix) #HPX_02298_RA_164.5_DEC_+38.7_STK.fits
		ps1_match = remove_ps1_source(ps1_cat, ns_cat)
		ns_cat = ns_cat[~ps1_match]
		#print len(ps1_match), ps1_match.sum()
		print "PS1 remove {0} stationary sources".format(ps1_match.sum())
		print "{0} non-stationary sources remain".format((~ps1_match).sum())
		fits.writeto('{}.ns'.format(fits_tab), ns_cat, header=sources.header)
	except ValueError:
		print 'ValueError: {0}'.format(fits_tab)
	except IOError:
		print 'IOError: {0}'.format(fits_tab)
	
def main():
	global RB_model
	global ps1_dir
	ps1_dir = '/mnt/jkavelaars/YTC_TEST/PS/'
	workdir = '/mnt/jkavelaars/YTC_TEST/'
	fits_list = filter(lambda x: x.endswith('fits'), os.listdir(workdir))
	fits_list = [workdir + i for i in fits_list]
	#fits_list = ['/sciproc/disk2/cfis/catalog_20170501/HPX_02559_RA_178.6_DEC_+35.7_cat.fits']
	training_fits = 'HPX_03011_RA_09.8_DEC_+30.0_cat.fits'
	RB_model = build_model(workdir + training_fits)
	map(ns_cat, fits_list)
	
if __name__ == '__main__':
	main()
	
	
	
	
	
	    
