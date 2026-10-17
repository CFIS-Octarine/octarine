from __future__ import absolute_import
import unittest
import mp_ephem
from astropy.io.fits import HDUList, PrimaryHDU

from daomop import downloader


class DownloaderTest(unittest.TestCase):
    """
    Simple test cases to ensure that downloader.py is returning the proper types for each of its methods.
    """
    def setUp(self):
        self.downloader = downloader.Downloader()
        self.obs = mp_ephem.ObsRecord.from_string(" c172f24340 * C2017 03 28.43594912 22 33.937+38 11 35.12         "
                                                  "23.5 r      568 O   2087679p34 c172f24340  Y  ------- ------- "
                                                  "0.20 4 23.54 0.10 % cand")

    def test_image_key(self):
        key = downloader.Downloader.image_key(self.obs)
        self.assertEqual(key, '2087679p34c172f24340')

    def test_get(self):
        get = self.downloader.get(self.obs)
        self.assertIsInstance(get, HDUList)

    def test_get_hdu(self):
        hdu_list = self.downloader.get_hdu(self.obs)
        for hdu in hdu_list:
            self.assertIsInstance(hdu, PrimaryHDU)
