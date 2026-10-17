from __future__ import absolute_import
import unittest
import mp_ephem
from daomop import storage


class FromFrameTest(unittest.TestCase):

    def setUp(self):
        self.filename = '../data/2002_VT130.0.ast'

    def test_return_from_string(self):
        tst = mp_ephem.ObsRecord.from_string(" 2002 VT130   C2017 03 24.22518 05 02 48.963+23 55 00.51                 "
                                             "    568 O   2086898p24 2002 VT130  Z   794.37 4347.80 0.20 4 ----- ---- %"
                                             " ")

        img = storage.FitsImage.from_frame(tst.comment.frame)

        self.assertIsInstance(img, storage.FitsImage)
        self.assertEquals(img.observation.dataset_name, "2086898")
        self.assertEquals(img.ccd, 24)
        self.assertEquals(img.uri, "vos:cfis/solar_system/dbimages/2086898/2086898p.fits.fz")

    def test_return_from_file(self):
        tst = mp_ephem.EphemerisReader()
        tst = tst.read(self.filename)
        img = storage.FitsImage.from_frame(tst[0].comment.frame)

        self.assertIsInstance(img, storage.FitsImage)
        self.assertEquals(img.observation.dataset_name, "2086898")
        self.assertEquals(img.ccd, 24)
        self.assertEquals(img.uri, "vos:cfis/solar_system/dbimages/2086898/2086898p.fits.fz")

    def test_incorrect_frame(self):
        short_dataset_name = mp_ephem.ObsRecord.from_string(" 2002 VT130   C2017 03 24.22518 05 02 48.963+23 55 00.51  "
                                                            "                   568 O   206898p24 2002 VT130  Z   "
                                                            "794.37 4347.80 0.20 4 ----- ---- % ")

        long_dataset_name = mp_ephem.ObsRecord.from_string(" 2002 VT130   C2017 03 24.22518 05 02 48.963+23 55 00.51   "
                                                           "                  568 O   20868298p24 2002 VT130  Z   "
                                                           "794.37 4347.80 0.20 4 ----- ---- % ")

        incorrect_ccd = mp_ephem.ObsRecord.from_string(" 2002 VT130   C2017 03 24.22518 05 02 48.963+23 55 00.51       "
                                                       "              568 O   2086898p54 2002 VT130  Z   794.37 4347.80"
                                                       " 0.20 4 ----- ---- % ")

        # dataset_name too short
        with self.assertRaises(ValueError):
            storage.FitsImage.from_frame(short_dataset_name.comment.frame)

        # dataset_name too long
        with self.assertRaises(ValueError):
            storage.FitsImage.from_frame(long_dataset_name.comment.frame)

        # ccd out of range
        with self.assertRaises(ValueError):
            storage.FitsImage.from_frame(incorrect_ccd.comment.frame)


if __name__ == '__main__':
    unittest.main()
