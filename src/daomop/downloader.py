import logging
import storage
from multiprocessing import Lock
from astropy.io import fits


class Downloader(object):
    def __init__(self):
        self._downloaded_images = dict()
        self.lock = Lock()
        self.locks = {}

    @staticmethod
    def image_key(obs_record):
        """
        Return the observation record's key.

        :param obs_record: the observation record whose key is returned
        :type obs_record: mp_ephem.ObsRecord
        :return: the key for the given observation record.
        """
        return obs_record.comment.frame + obs_record.provisional_name

    def get(self, obs_record):
        """
        Returns the HDUList of a given observation record.
        Uses locks on threads for multiprocessing to retrieve multiple HDU's in parallel.

        :param obs_record: the ObsRecord for which the image is to be retrieved
        :type obs_record: mp_ephem.ObsRecord
        :return: astropy.io.fits.hdu.image.HDUList
        """
        with self.lock:
            if self.image_key(obs_record) not in self.locks:
                self.locks[self.image_key(obs_record)] = Lock()

        with self.locks[self.image_key(obs_record)]:
            if self.image_key(obs_record) not in self._downloaded_images:
                self._downloaded_images[self.image_key(obs_record)] = self.get_hdu(obs_record)

        return self._downloaded_images[self.image_key(obs_record)]

    @staticmethod
    def put(artifact):
        """
        Put the artifact to its database.

        :param artifact: the file object to store to VOSpace
         :type artifact: storage.Artifact
        """
        # Lock off other network actions while upload in progress.
        artifact.put()

    def get_hdu(self, obs_record):
        """
        Retrieve a fits image associated with a given obs_record and return the list of HDUs off the associated cutout.

        :param obs_record: the ObsRecord for which the image is to be retrieved
        :type obs_record: mp_ephem.ObsRecord
        :return: [fits.ImageHDU, ... ]
        """
        logging.debug("Retrieving {}".format(self.image_key(obs_record)))
        try:
            image = storage.FitsImage.from_frame(obs_record.comment.frame)
            hdu_list = image.ra_dec_cutout(obs_record.coordinate, trim_to_datasec=True)
        except Exception as ex:
            logging.debug(ex)
            raise ex

        if not isinstance(hdu_list, fits.HDUList):
            return None

        logging.debug("Got {}".format(self.image_key(obs_record)))
        return hdu_list
