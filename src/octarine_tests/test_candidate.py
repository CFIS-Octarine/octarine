from __future__ import absolute_import
import unittest
from copy import deepcopy

from daomop import candidate


class CandidateTest(unittest.TestCase):
    """
    Testing CandidateSet object from candidate.py. Test revolve around the next() method which returns the next
     set of candidates.
    """
    def setUp(self):
        # assuming that vos:/cfis/solar_system/catalogs/17AQ06 still exists
        self.candidates = candidate.CandidateSet(2434, catalog_dir='17AQ06')

    def test_next(self):

        candidates = deepcopy(self.candidates)

        self.assertEqual(candidates.target.current_observation, -1)
        self.assertEqual(len(candidates.target.observation_sets), 0)
        self.assertEqual(candidates.catalog.current_date, 0)

        candidates.next()

        self.assertEqual(candidates.target.current_observation, 0)
        self.assertEqual(len(candidates.target.observation_sets), 29)
        self.assertEqual(candidates.catalog.current_date, 1)

        self.assertEqual(candidates.catalog.catalog.dataset_name, 'HPX_02434_RA_185.6_DEC_+37.2')
        self.assertEqual(candidates.catalog.catalog.filename, 'HPX_02434_RA_185.6_DEC_+37.2_bk.json')

        provisional_name = candidates.target.provisional_name.split(' ')[0]
        self.assertEqual(provisional_name, 'c172f24340')

    def test_provisional(self):
        candid = deepcopy(self.candidates.next())
        self.assertEqual(candid[0].provisional_name, 'c172f24340')
