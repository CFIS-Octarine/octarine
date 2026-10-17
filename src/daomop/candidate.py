import sys

from mp_ephem import ObsRecord
from astropy.time import Time
import logging
from . import storage

_LETTERS = 'abcdefghijklmnopqrstuvwxyz'
_DIGITS = '0123456789'


def provisional(mjd, hpx, count):
    """Compute a provisional name from the date and HPX information.
    """
    discovery_date = Time(int(mjd), format='mjd').yday.split(':')
    discovery_day = int(discovery_date[1])
    discovery_year = int(discovery_date[0])

    yr = discovery_year - 2000
    p1 = discovery_day//36
    p2 = discovery_day - p1*36
    dy = (_DIGITS + _LETTERS)[p1] + (_DIGITS + _LETTERS)[p2]
    p1 = count
    s2 = ""
    logging.info("Loading candidate : {}".format(count))
    while p1 > 35:
        p2 = p1//36
        s2 += (_DIGITS + _LETTERS)[p1 - p2 * 36]
        p1 = p2
    s1 = (_DIGITS + _LETTERS)[p1]
    return "c{:02}{:2}{:04}{:1}{:1}".format(yr, dy, hpx, s1, s2)


class ObservationSet(object):
    def __init__(self, provisional_name, record):
        self.provisional_name = provisional_name
        self.record = record
        self.current_observation = -1

    def __iter__(self):
        return self

    def next(self):
        """

        :return: An Observation Record
         :rtype: ObsRecord
        """
        if not self.current_observation + 1 < len(self.record['mag']):
            raise StopIteration
        self.current_observation += 1
        return ObsRecord(provisional_name=self.provisional_name,
                         discovery=True,
                         note1=None,
                         note2='C',
                         date=Time(self.record['mjd'][self.current_observation], format='mjd'),
                         ra=self.record['ra'][self.current_observation],
                         dec=self.record['dec'][self.current_observation],
                         mag=self.record['mag'][self.current_observation],
                         band=self.record['filterid'][self.current_observation],
                         observatory_code=568,
                         mag_err=self.record['magerr'][self.current_observation],
                         comment='cand',
                         xpos=None,
                         ypos=None,
                         frame=self.record['fitsname'][self.current_observation],
                         plate_uncertainty=None,
                         astrometric_level=4
                         )


class Target(object):
    def __init__(self, hpx, mjdate, record):
        self.hpx = hpx
        self.mjdate = mjdate
        self.record = record
        self.current_observation = -1

    def __iter__(self):
        return self

    @property
    def observation_sets(self):
        return self.record.keys()

    @property
    def provisional_name(self):
        return provisional(self.mjdate, self.hpx, self.current_observation)

    def next(self):
        """

        :rtype: ObservationSet
        """
        self.current_observation += 1
        if not self.current_observation < len(self.observation_sets):
            raise StopIteration
        return ObservationSet(self.provisional_name,
                              self.record[self.observation_sets[self.current_observation]])

    def previous(self):
        if self.current_observation == 0:
            return None
        self.current_observation -= 1
        return ObservationSet(self.provisional_name,
                              self.record[self.observation_sets[self.current_observation]])


class Catalog(object):
    def __init__(self, pixel, catalog_dir=None):
        self.catalog = storage.JSONCatalog(pixel, catalog_dir=catalog_dir)
        self.current_date = -1
        self.current_target = -1

    def __iter__(self):
        return self

    @property
    def mjdates(self):
        return self.catalog.json.keys()

    def next(self):
        """

        :rtype: Target
        """
        self.current_date += 1
        if not self.current_date < len(self.mjdates):
            raise StopIteration
        mjdate = self.mjdates[self.current_date]
        return Target(self.catalog.pixel, mjdate, self.catalog.json[mjdate])

    def previous(self):
        if self.current_date == 0:
            return
        self.current_date -= 1
        mjdate = self.mjdates[self.current_date]
        return Target(self.catalog.pixel, mjdate, self.catalog.json[mjdate])


class CandidateSet(object):

    def __init__(self, pixel, catalog_dir=None):
        self.catalog = Catalog(pixel, catalog_dir=catalog_dir)
        self.target = self.catalog.next()

    def __iter__(self):
        return self

    def next(self):
        """

        :rtype: list(ObsRecord)
        """
        try:
            observation_set = self.target.next()
            obs = []
            for ob in observation_set:
                obs.append(ob)
            return obs
        except StopIteration:
            self.target = self.catalog.next()
            return self.next()

    def previous(self):
        try:
            observation_set = self.target.previous()
            if observation_set is None:
                return None
            obs = []
            # noinspection PyTypeChecker
            for ob in observation_set:
                obs.append(ob)
            return obs
        except StopIteration:
            self.target = self.catalog.previous()
            return self.previous()

if __name__ == "__main__":
    print sys.argv
    for target in Catalog(int(sys.argv[1])):
        for obs_set in target:
            for observation in obs_set:
                print(observation)
            print("")
