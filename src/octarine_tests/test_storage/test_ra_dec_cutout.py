from __future__ import absolute_import
import unittest
import mp_ephem
import hashlib
from daomop import storage
from astropy.io import fits

TEST_INTEGER = 10
TEST_FLOAT = 1.0/60.0


class FitsImageGetTest(unittest.TestCase):
    """
    Tests that proper errors are thrown and return types are correct in ra_dec_cutout method.
    Can take 10+ seconds to complete a full test due to file downloads from VOSpace.
    """
    def setUp(self):
        self.obs = mp_ephem.ObsRecord.from_string(" 2002 VT130   C2017 03 24.22518 05 02 48.963+23 55 00.51            "
                                                  "         568 O   2086898p24 2002 VT130  Z   794.37 4347.80 0.20 4 --"
                                                  "--- ---- % ")

        storage.DBIMAGES = "vos:jkavelaars/TNORecon/dbimages"

    def test_TypeError(self):
        fit = storage.FitsImage.from_frame(self.obs.comment.frame)

        with self.assertRaises(TypeError):
            fit.ra_dec_cutout(fit.ccd)

        with self.assertRaises(TypeError):
            fit.ra_dec_cutout(self.obs.coordinate, TEST_INTEGER)

        with self.assertRaises(TypeError):
            fit.ra_dec_cutout(self.obs.coordinate, self.obs)

    def test_ReturnType(self):
        self.assertIsInstance(storage.FitsImage.from_frame(self.obs.comment.frame).ra_dec_cutout(self.obs.coordinate),
                              fits.hdu.HDUList)

        self.assertIsInstance(storage.FitsImage.from_frame(self.obs.comment.frame).ra_dec_cutout(self.obs.coordinate,
                                                                                                 1.0/60.0),
                              fits.hdu.HDUList)

    def test_MD5Sum(self):
        m = hashlib.md5()
        m.update(storage.FitsImage.from_frame(self.obs.comment.frame).ra_dec_cutout(self.obs.coordinate)[1].data)

        # written on 2017/05/31, assumes cutout() and other associated methods work properly
        self.assertEqual(m.hexdigest(), "69e8f3f49f000ac78ea83db63468095c")


if __name__ == '__main__':
    unittest.main()
