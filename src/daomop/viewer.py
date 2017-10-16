import os
import re
import logging
from multiprocessing.dummy import Pool, Lock
from math import atan2, degrees
from multiprocessing.pool import ApplyResult

from astropy.time import Time
from cadcutils.exceptions import NotFoundException
from mp_ephem import ObsRecord

from . import candidate
from . import downloader
from . import storage
from ginga import AstroImage
from ginga.web.pgw import ipg, Widgets, Viewers
from ginga.misc import log
from astropy.wcs import WCS

logging.basicConfig(level=logging.INFO, format="%(module)s.%(funcName)s:%(lineno)s %(message)s")
DISPLAY_KEYWORDS = ['EXPNUM', 'DATE-OBS', 'UTC-OBS', 'EXPTIME', 'FILTER']
LEGEND = 'Keyboard Shortcuts: \n' \
         'f: image backwards \n' \
         'g: image forwards \n' \
         'q: pan mode \n(click and drag on canvas)\n' \
         't: contrast mode \n(right click on canvas after pressing "t" to reset contrast)\n' \
         'esc: reset keyboard mode\n'
ACCEPTED_DIRECTORY = 'accepted'
PROCESSES = 5


class ConsoleBoxStream(object):
    """
    A class that writes to a console box as a stream.
    """
    def __init__(self, console_box):
        self.console_box = console_box

    def write(self, content):
        """

        :param content: content to write to the console stream.
        :return:
        """
        self.console_box.append_text(content)

    def flush(self):
        pass


class ValidateGui(ipg.EnhancedCanvasView):

    def __init__(self, logger, window, bindings=None):
        """

        :param logger: a logger object to send messages to
        :type logger: logging.Logger
        :param window: The main window of the application
        :param bindings: Any bindings previously set on this window.
        """
        super(ValidateGui, self).__init__(logger=logger, bindings=bindings)

        self.console_box = Widgets.TextArea(editable=False)

        self.downloader = downloader.Downloader()
        self.pool = Pool(processes=PROCESSES)
        self.lock = Lock()
        self.image_list = {}
        self.astro_images = {}

        self.logger = logger
        console_handler = logging.StreamHandler(stream=ConsoleBoxStream(self.console_box))
        self.logger.addHandler(console_handler)
        self.top = window

        self.enable_autocuts('on')
        self.set_autocut_params('zscale')

        # creating drawing canvas; initializing polygon types
        self.canvas = self.add_canvas()
        self.circle = self.canvas.get_draw_class('circle')

        # creating key-press event handling
        self.canvas.add_callback('key-press', self._key_press, 'key', self)

        # remove callbacks for clicking on the canvas (which is over the viewer)
        self.canvas.delete_callback('cursor-down')
        self.canvas.delete_callback('cursor-up')

        self.obs_number = 0
        self.candidate = None
        self.candidates = None
        self.zoom = None
        self._center = None
        self.healpix = None
        self.storage_list = None
        self.override = None
        self.qrun_id = None
        self.length_check = False

        # GUI elements
        self.pixel_base = 1.0
        self.readout = Widgets.Label("")
        self.header_box = Widgets.TextArea(editable=False)
        self.accept = Widgets.Button("Accept")
        self.reject = Widgets.Button("Reject")
        self.next_set = Widgets.Button("Next Set >")
        self.previous_set = Widgets.Button("< Previous Set")
        self.load_json = Widgets.Button("Load")
        self.clear_button = Widgets.Button("Clear")
        self.yes_button = Widgets.Button("Yes")
        self.no_button = Widgets.Button("No")
        self.reload_button = Widgets.Button("Reload")
        self.warning = Widgets.Label("In case you try to reject a previously accepted candidate: ")

        self.legend = Widgets.TextArea(wrap=True)
        self.legend.set_text(LEGEND)
        self.build_gui(self.top)
        self.comparison_images = {}
        self.null_observation = {}
        self.next_image = None

    def build_gui(self, container):
        """
        Building the GUI to be displayed in an HTML5 canvas.
        Tested and working in Mozilla Firefox and Google Chrome web browsers.

        :param container: ginga.web.pgw.Widgets.TopLevel object
        """
        bindings = self.get_bindings()
        bindings.enable_all(True)

        # keyboard mode indicator, upper right corner
        self.show_mode_indicator(True, corner='ur')

        viewer_vbox = Widgets.VBox()  # box containing the viewer
        viewer_vbox.set_border_width(2)
        viewer_vbox.set_spacing(1)
        viewer_widget = Viewers.GingaViewerWidget(viewer=self)
        viewer_vbox.add_widget(viewer_widget, stretch=1)
        viewer_vbox.add_widget(self.readout, stretch=0)  # text directly below the viewer for coordinate display

        self.set_callback('cursor-changed', self.motion_cb)

        healpix_set = Widgets.TextEntrySet()
        healpix_set.add_callback('activated', lambda x: self.set_healpix(event=x))
        healpix_set.set_length(6)

        candidate_override = Widgets.TextEntrySet()
        candidate_override.add_callback('activated', lambda x: self.override_set(event=x))
        candidate_override.set_length(10)

        astfile = Widgets.TextEntry(editable=True)
        astfile.add_callback('activated', lambda x: self.load_astfile(event=x))

        catalog = Widgets.TextEntrySet(text='17AQ03')
        catalog.add_callback('activated', lambda x: self.set_qrun_id(x))
        catalog.set_length(5)

        self.accept.add_callback('activated', lambda x: self.accept_reject())
        self.reject.add_callback('activated', lambda x: self.accept_reject(rejected=True))
        self.load_json.add_callback('activated', lambda x: self.load_candidates())
        self.next_set.add_callback('activated', lambda x: self.next())
        self.previous_set.add_callback('activated', lambda x: self.previous())
        self.clear_button.add_callback('activated', lambda x: self.clear_viewer())
        self.reload_button.add_callback('activated', lambda x: self.reload_candidates())

        # accept/reject/next buttons
        buttons_hbox = Widgets.HBox()
        buttons_hbox.add_widget(self.previous_set)
        buttons_hbox.add_widget(self.accept)
        buttons_hbox.add_widget(self.reject)
        buttons_hbox.add_widget(self.next_set)
        buttons_hbox.add_widget(self.load_json)
        buttons_hbox.add_widget(self.reload_button)
        buttons_hbox.add_widget(self.clear_button)
        self.load_json.set_enabled(False)
        buttons_hbox.set_spacing(3)

        # catalog directory text box
        catalog_box = Widgets.HBox()
        catalog_label = Widgets.Label(text="Set QRUNID:", style='color:red')
        catalog_box.add_widget(catalog_label)
        catalog_box.add_widget(catalog)
        catalog_box.set_margins(15, 0, 10, 0)  # top, right, bottom, left

        candidates_hbox = Widgets.HBox()
        candidate_label = Widgets.Label(text="(Optional) Enter candidate set (HEALPIX): ")
        candidates_hbox.add_widget(candidate_label)
        candidates_hbox.add_widget(healpix_set)
        candidates_hbox.set_margins(15, 0, 15, 0)  # top, right, bottom, left

        override_hbox = Widgets.HBox()
        override_label = Widgets.Label(text="(Optional) Override provisional name: ")
        override_hbox.add_widget(override_label)
        override_hbox.add_widget(candidate_override)
        override_hbox.set_margins(0, 0, 15, 0)  # top, right, bottom, left

        astfile_hbox = Widgets.HBox()
        astfile_hbox_label = Widgets.Label(text="Paste AST file here:")
        astfile_hbox.add_widget(astfile_hbox_label)
        astfile_hbox.add_widget(astfile)

        # button and text entry vbox
        buttons_vbox = Widgets.VBox()
        buttons_vbox.add_widget(buttons_hbox)
        buttons_vbox.add_widget(catalog_box)
        buttons_vbox.add_widget(candidates_hbox)
        buttons_vbox.add_widget(override_hbox)
        buttons_vbox.add_widget(astfile_hbox)

        viewer_vbox.add_widget(buttons_vbox)  # add buttons below the viewer

        viewer_header_hbox = Widgets.HBox()  # box containing the viewer/buttons and rightmost text area
        viewer_header_hbox.add_widget(viewer_vbox)
        viewer_header_hbox.add_widget(Widgets.Label(''))
        hbox = Widgets.HBox()
        hbox.add_widget(self.header_box)
        hbox.add_widget(self.legend)
        viewer_header_hbox.add_widget(hbox)

        full_vbox = Widgets.VBox()  # vbox container for all elements
        full_vbox.add_widget(viewer_header_hbox)

        full_vbox.add_widget(self.console_box)
        self.console_box.set_text('Logging output:\n')
        self.header_box.set_text("Header:")

        container.set_widget(full_vbox)
        container.set_widget(self.warning)
        container.set_widget(self.yes_button)
        container.set_widget(self.no_button)
        self.yes_button.set_enabled(False)
        self.no_button.set_enabled(False)
        self.buttons_off()

    def next(self):
        """
        Load the next set of images into the viewer
        """
        if self.candidates is not None:
            # noinspection PyBroadException
            try:
                self.buttons_off()
                self.obs_number = 0
                self.next_image = None
                self.clear_candidate_images()
                self.candidate = self.candidates.next()

                # finding next candidate to load depending on which .ast files are written
                while self.ast_exists():
                    self.candidate = self.candidates.next()

                self.logger.info("Loading {}...".format(self.candidate[0].provisional_name))
                self.load()
                self.buttons_on()

            except Exception as ex:
                self.logger.info('Loading next candidate set failed.')
                if isinstance(ex, StopIteration):
                    self.logger.info('StopIteration error: End of candidate set.')
                    self.logger.info('Hit "Load" button to move onto the next set.')
                    self.previous_set.set_enabled(True)
                    self.load_json.set_enabled(True)

    def ast_exists(self):
        """
        Checks if candidate has already been examined with a file written on VOSpace.
        length_check is necessary because it means the sub directory exists, if it doesn't an error will be
         thrown when looking in the directory list.

        :return True is the .ast files exists and there is no viewing override, False otherwise
        """
        if self.length_check and self.candidate[0].provisional_name + '.ast' in storage.listdir(
                os.path.join(os.path.dirname(storage.DBIMAGES), storage.CATALOG, self.qrun_id,
                             self.candidates.catalog.catalog.dataset_name), force=True):

            if self.override == self.candidate[0].provisional_name:
                self.logger.info("Candidate {} being overridden for viewing."
                                 .format(self.candidate[0].provisional_name))

            else:
                self.logger.info("Candidate {} has been investigated.".format(self.candidate[0].provisional_name))
                return True

        return False

    def previous(self):
        """
        Load the previous set of images into the viewer
        """
        if self.candidates is not None:
            self.buttons_off()
            self.obs_number = 0
            self.candidate = self.candidates.previous()
            if self.candidate is not None:
                self.load()
            self.buttons_on()

    def accept_reject(self, rejected=False):
        """
        Accept or reject current observation depending on button press. Write to file and load next set into the viewer

        :param rejected: whether the candidate set has been accepted or rejected
        """
        self.logger.info("Rejected.") if rejected else self.logger.info("Accepted.")
        self.buttons_off()
        if self.candidates is not None:
            self.write_record(rejected=rejected)
            self.next()

    def set_qrun_id(self, qrun_id):
        """
        :param qrun_id: QRUNID in a header file
        """
        if hasattr(qrun_id, 'text'):
            self.qrun_id = str(qrun_id.text).strip(' ')
            self.storage_list = storage.listdir(os.path.join(os.path.dirname(storage.DBIMAGES),
                                                             storage.CATALOG,
                                                             self.qrun_id), force=True)
            self.load_json.set_enabled(True)
            self.logger.info("QRUNID set to {}.".format(self.qrun_id))

    def lookup(self):
        """
        Determines which healpix values is to be examined next. The healpix value is eventually used when creating the
        candidate set.

        :return healpix value for the candidate files; 0 if no candidate files have been found
        """
        count = 0
        self.length_check = False
        for filename in self.storage_list:

            # ex: filename = HPX_00887_RA_203.6_DEC_+58.9_bk.json,
            #     filename[:-len(storage.MOVING_TARGET_VERSION)] = HPX_00887_RA_203.6_DEC_+58.9
            # sub_directory will be the directory where a candidate's .ast files are written
            sub_directory = filename[:-len(storage.MOVING_TARGET_VERSION + '.json')]
            count += 1

            # if the file extension is in the filename, then it is a file containing candidate information
            if storage.MOVING_TARGET_VERSION in filename:
                x = re.match('(?P<hpx>HPX_)(?P<healpix>\d{5})(?P<leftover>_.*)', filename)

                if self.healpix is not None and int(x.group('healpix')) < self.healpix:
                    continue  # skipping over json files until the specified catalog has been reached

                # if the sub directory exists, we will have to check that all the candidates have been investigated
                elif sub_directory in self.storage_list:
                    self.length_check = True

                # TODO: go back to server for storage_list in case two people are actively writing from unique servers
                # cutting down the storage list for further iterating
                self.storage_list = self.storage_list[count:]
                return int(x.group('healpix'))
        return 0

    def set_healpix(self, event):
        """
        Sets the healpix for the current Candidate set.

        :param event: healpix value
        """
        if hasattr(event, 'text'):
            self.healpix = int(event.text)
            self.logger.info("Set healpix as {}".format(self.healpix))
            if self.qrun_id is not None:
                self.load_json.set_enabled(True)

    def load_astfile(self, event):
        self.candidate = []
        for line in event.text.split('\n'):
            # noinspection PyBroadException
            try:
                obs_record = ObsRecord.from_string(line)
                if obs_record is not None:
                    self.candidate.append(obs_record)
            except:
                logging.warning("Failed to parse line >{}<".format(line))
                return
        self.logger.info("Accepted AST file.")
        self.candidates = [self.candidate]
        self.next_set.set_enabled(False)
        self.previous_set.set_enabled(False)
        self._download_obs_records(self.candidate)
        self.load(0)

    def override_set(self, event):
        """
        Look at the cutout even if it has already been investigated. Primarily used for double checking
         accepted candidates.
        """
        if hasattr(event, 'text'):
            self.override = str(event.text).strip(' ')
            self.logger.info("Will override {}.".format(self.override))

    def load_candidates(self, healpix=None):
        """
        Initial candidates loaded into the viewer. Starts up a thread pool to download images simultaneously.

        :param healpix: Catalogue number containing dataset
        """
        if healpix is None:
            self.healpix = self.lookup()

        self.buttons_off()

        while self.healpix != 0 and self.set_examined():
            self.healpix = self.lookup()
            if self.healpix == 0:  # base case (when there are no more open candidate sets in the VOSpace directory)
                self.logger.info("No more candidate sets for this QRUNID.")
                raise StopIteration

        self.logger.warning("Launching image prefetching. Please be patient.")

        with self.lock:
            for obs_records in self.candidates:
                self._download_obs_records(obs_records)

        self.candidates = candidate.CandidateSet(self.healpix, catalog_dir=self.qrun_id)
        self.candidate = None  # reset on candidate to clear it of any leftover from previous sets
        self.load()

    def _download_obs_records(self, record):
        """
        Download the observations associated with the current self.candidate set of obsRecords.
        :return:
        """
        previous_record = None
        previous_offset = 2 * storage.CUTOUT_RADIUS
        offset = previous_offset
        for obs_record in record:
            assert isinstance(obs_record, ObsRecord)
            key = self.downloader.image_key(obs_record)
            if key not in self.image_list:
                self.image_list[key] = self.pool.apply_async(self.downloader.get, (obs_record,))

            # Check if we should load a comparison for the previous image.
            if previous_record is not None:
                offset = obs_record.coordinate.separation(previous_record.coordinate)
                if offset > storage.CUTOUT_RADIUS and previous_offset > storage.CUTOUT_RADIUS:
                    # Insert a blank image in the list
                    previous_key = self.downloader.image_key(previous_record)
                    comparison = storage.get_comparison_image(previous_record.coordinate,
                                                              previous_record.date.mjd)
                    frame = "{}{}".format(comparison[0]['observationID'], 'p00')
                    comparison_obs_record = ObsRecord(null_observation=True,
                                                      provisional_name=previous_record.provisional_name,
                                                      date=Time(comparison[0]['mjdate'], format='mjd',
                                                                precision=5).mpc,
                                                      ra=previous_record.coordinate.ra.degree,
                                                      dec=previous_record.coordinate.dec.degree,
                                                      frame=frame,
                                                      comment=previous_key)
                    key = self.downloader.image_key(comparison_obs_record)
                    self.null_observation[key] = comparison_obs_record
                    self.comparison_images[previous_key] = key
                    if key not in self.image_list:
                        self.image_list[key] = self.pool.apply_async(self.downloader.get,
                                                                     (comparison_obs_record,))

            previous_record = obs_record
            previous_offset = offset
        # Check if the offset between the last record and the one just before it was large.
        if previous_offset > storage.CUTOUT_RADIUS and previous_record is not None:
            previous_key = self.downloader.image_key(previous_record)
            comparison = storage.get_comparison_image(previous_record.coordinate,
                                                      previous_record.date.mjd)
            frame = "{}{}".format(comparison[0]['observationID'], 'p00')
            comparison_obs_record = ObsRecord(null_observation=True,
                                              provisional_name=previous_record.provisional_name,
                                              date=Time(comparison[0]['mjdate'], format='mjd',
                                                        precision=5).mpc,
                                              ra=previous_record.coordinate.ra.degree,
                                              dec=previous_record.coordinate.dec.degree,
                                              frame=frame,
                                              comment=previous_key)
            key = self.downloader.image_key(comparison_obs_record)
            self.null_observation[key] = comparison_obs_record
            self.comparison_images[previous_key] = key
            if key not in self.image_list:
                self.image_list[key] = self.pool.apply_async(self.downloader.get,
                                                             (comparison_obs_record,))

    def set_examined(self):
        """
        Checks if the current json file has been fully examined or not

        :return True if the directory is fully examined and there's no override, False if it has not been examined.
        """
        self.logger.info("Accepted candidate entry: {}".format(self.healpix))
        try:
            self.candidates = candidate.CandidateSet(self.healpix, catalog_dir=self.qrun_id)
            if self.length_check:
                sub_directory = storage.listdir(os.path.join(os.path.dirname(storage.DBIMAGES),
                                                             storage.CATALOG,
                                                             self.qrun_id,
                                                             self.candidates.catalog.catalog.dataset_name), force=True)
                if self.override is not None:
                    filename = self.override+'.ast'
                    if filename in sub_directory:
                        self.logger.info("Overriding {}.".format(filename))
                        return False
                else:
                    count = 0
                    # counting the total amount of candidates that are in self.candidates
                    for _ in self.candidates:
                        count += 1

                    # re-set self.candidates since the for loop removes all its candidates in a dequeuing fashion
                    self.candidates = candidate.CandidateSet(self.healpix, catalog_dir=self.qrun_id)

                    # the amount of files in the accompanying subdirectory for the .json candidate file
                    directory_length = len(sub_directory)
                    if count == directory_length:
                        self.logger.info("Candidate set {} fully examined.".format(self.healpix))
                        return True

                    elif count > directory_length:
                        self.logger.info("Candidate set {} not fully examined.".format(self.healpix))
                        return False

                    else:
                        self.logger.error("Value error: count {} or directory_length {} is out of range."
                                          .format(count, directory_length))
                        raise ValueError

            return False  # no length check, therefor no directory has been created and this set isn't examined

        except Exception as ex:
            self.logger.info("Failed to load candidates: {}".format(str(ex)))
            if isinstance(ex, StopIteration):
                self.logger.info('StopIteration error. Candidate set might be empty.')
                return True  # continue with iteration
            else:
                raise ex

    def reload_candidates(self):
        """
        Performs a hard reload on all images for the case of loading errors.
        Closes current worker pool and reopens a new one.
        """
        if self.healpix is not None:
            self.logger.info('Reloading all candidates...')
            self.pool.terminate()
            self.pool = Pool(processes=PROCESSES)
            self.buttons_on()
            self.set_qrun_id(self.qrun_id)
            self.load_candidates(self.healpix)
            self.next()

    def load(self, obs_number=0):
        """
        With the viewing window already created, Creates a FitsImage object and loads its cutout into the window and
         displays select header values (see: DISPLAY_KEYWORDS).
        Define the center of the first image to be the reference point for aligning the other two images in the set.

        :param obs_number: index of which line in the file gets loaded/displayed in the viewer
        """
        self._center = None
        self.obs_number = obs_number
        self._load()
        self._center = WCS(self.header).all_pix2world(self.get_data_size()[0] / 2,
                                                      self.get_data_size()[1] / 2, 0)

    def _load(self):
        """
        Loads an image into the viewer, applying appropriate transformations for proper display.
        Checks if an HDU has been loaded already and retrieves if needed and then displays that HDU.
        Uses multiprocessing techniques for simultaneous downloads and dictionaries to keep track of which images
         have been already loaded for faster image switching.
        """
        # load the image if not already available, for now we'll put this in here.
        if self.candidates is None:
            self.logger.info("No candidates loaded.")
            return

        # loads first candidate
        if self.candidate is None:
            self.next()
            return

        key = self.key
        while True:
            # noinspection PyBroadException
            try:
                if key not in self.astro_images:
                    # TODO: MEF
                    image = AstroImage.AstroImage(logger=self.logger)
                    image.load_hdu(self.loaded_hdu)
                    self.astro_images[key] = image

                self.set_image(self.astro_images[key])

                if self.zoom is not None:
                    self.zoom_to(self.zoom)
                self._rotate()
                if self.center is not None:
                    self._align()

                # the image cutout is considered the first object on the canvas, this deletes everything over top of it
                self.canvas.delete_objects(self.canvas.get_objects()[1:])
                if key not in self.null_observation:
                    self.mark_aperture()

                self.header_box.set_text("Header:" + self.info)
                self.logger.info("Loaded: {}".format(self.candidate[self.obs_number].comment.frame))
                break

            except Exception as ex:
                self.logger.info(str(ex))
                self.logger.info("Skipping candidate {} due to load failure."
                                 .format(self.candidate[0].provisional_name))
                self.next()
                break

    def mark_aperture(self):
        """
        Draws a red circle on the drawing canvas in the viewing window around the celestial object detected.
        """
        ra = self.candidate[self.obs_number].coordinate.ra
        dec = self.candidate[self.obs_number].coordinate.dec
        x, y = WCS(self.header).all_world2pix(ra, dec, 0)
        self.canvas.add(self.circle(x, y, radius=10, color='red'))

    def write_record(self, rejected=False):
        """
        Writing observation lines to a new file.

        :param rejected: Whether or not the candidate set contains a valid celestial object
        :type rejected: bool
        """
        try:
            catalog_dir = os.path.join(storage.CATALOG,
                                       self.qrun_id,
                                       self.candidates.catalog.catalog.dataset_name)

            art = storage.ASTRecord(self.candidate[0].provisional_name,
                                    version='',
                                    catalog_dir=catalog_dir)

            with open(art.filename, 'w+') as fobj:
                for ob in self.candidate:
                    if rejected:
                        ob.null_observation = True
                    fobj.write(ob.to_string() + '\n')

            self.logger.info("Queuing job to write file to VOSpace.")
            with self.lock:
                try:
                    if rejected:
                        self.remove_check(art)
                    elif not rejected:
                        self.pool.apply_async(self.accepted_list, (art,))
                        self.logger.info(
                            "Done Queuing {} for VOSpace write.".format(self.candidate[0].provisional_name + ".ast"))

                except Exception as ex:
                    self.logger.info("Failed to write file {}: {}".format(
                        self.candidate[0].provisional_name, str(ex)))

        except IOError as ex:
            self.logger.info("Unable to write to file.")
            self.logger.info(str(ex))
            raise ex

    def remove_check(self, art, ext='.ast'):
        """
        Checks a file's existence in its /accepted/ VOSpace directory. If the uri can't be found, there must not be
         a file in it's accepted directory, so a standard null observation file is uploaded.
        Prompts user to make a decision if the file already exists.

        :param art: artifact who's being checked for existence
        :param ext: file type
        """
        # noinspection PyBroadException
        try:
            accepted_uri = os.path.join(os.path.join(os.path.dirname(storage.DBIMAGES), storage.CATALOG),
                                        self.header['QRUNID'], ACCEPTED_DIRECTORY, art.observation.dataset_name + ext)
            if storage.exists(accepted_uri):
                self.yes_button.add_callback('activated', lambda x: self.move_accepted(accepted_uri, art))
                self.no_button.add_callback('activated', lambda x: self.warning_label_reset())
                self.yes_button.set_enabled(True)
                self.no_button.set_enabled(True)

                self.logger.warning("File already accepted.")
                self.warning.set_text("FILE {} HAS ALREADY BEEN ACCEPTED, ARE YOU SURE YOU WANT TO REJECT IT?"
                                      .format(art.observation.dataset_name))

                self.warning.set_color(fg='white', bg='red')

        except Exception as ex:
            if isinstance(ex, NotFoundException):
                self.write_rejected(art)
            else:
                self.logger.error(str(ex))
                raise ex

    def move_accepted(self, accepted_uri, art):
        """
        Deletes the file at the uri and queue's a thread to write the file in a new destination as a rejected
        observation. Disables buttons and resets label.

        :param accepted_uri: uri of the accepted file
        :param art: artifact object for the record being examined
        """
        storage.delete(accepted_uri)
        self.logger.info("Deleted {}".format(accepted_uri))
        self.write_rejected(art)
        self.warning_label_reset()

    def warning_label_reset(self):
        """
        Method that serves as a callback destination. Disables yes/no buttons and resets label text.
        """
        self.yes_button.set_enabled(False)
        self.no_button.set_enabled(False)
        self.warning.set_text("In case you try to reject a previously accepted candidate: ")
        self.warning.set_color(fg='black', bg='white')

    def write_rejected(self, art):
        """
        Start a thread to write the rejected artifact to its uri

        :param art: Artifact object
        """
        self.pool.apply_async(self.downloader.put, (art,))
        self.logger.info("Done Queuing {} for VOSpace write {}".format(self.candidate[0].provisional_name + ".ast",
                                                                       art.uri))

    def accepted_list(self, art, ext='.ast'):
        """
        Places accepted .ast file in an accepted folder in its QRUNID section on VOSpace

        :param art: Artifact object containing the proper file name
        :param ext: file extension
        """
        # 'vos:cfis/solar_system/dbimages/catalogs/<QRUNID>/accepted/<dataset_name>.ast
        # Since this just uploads an unintuitive name in the directory, perhaps the path could be changed to
        #  ../accepted/<healpix>/<dataset_name>.ast
        destination = os.path.join(os.path.join(os.path.dirname(storage.DBIMAGES), storage.CATALOG),
                                   self.header['QRUNID'], ACCEPTED_DIRECTORY, art.observation.dataset_name + ext)
        try:
            storage.make_path(destination)
            storage.copy(art.filename, destination)
        except Exception as ex:
            self.logger.info("Failed writing to accepted directory for {}: {}"
                             .format(art.observation.dataset_name, str(ex)))
            raise ex

    def _rotate(self):
        """
        Rotates the current viewer image to be oriented North up East left. This is done by taking outward vectors from
         the origin and using their WCS values to determine the original orientation of the image. Images are then
         flipped/rotated accordingly to be North up East left.
        """
        wcs = WCS(self.header)
        self.transform(False, False, False)
        x = wcs.all_pix2world([[0, 0], [1, 1], [1, 0]], 0)
        ra1 = x[0][0]
        ra2 = x[1][0]
        ra3 = x[2][0]
        dec1 = x[0][1]
        dec2 = x[1][1]
        dec3 = x[2][1]

        delta_x = ra2 - ra1
        delta_y = dec2 - dec1

        flip_x = 1
        flip_y = 1
        if not delta_x < 0:
            flip_x = -1
            if not delta_y > 0:
                flip_y = -1
                self.transform(True, True, False)  # def transform(self, flip_x, flip_y, swap_xy):
            else:
                self.transform(True, False, False)
        elif not delta_y > 0:
            flip_y = -1
            self.transform(False, True, False)

        delta_delta = (dec3 - dec1) * flip_y
        delta_ra = (ra1 - ra3) * flip_x

        theta = degrees(atan2(delta_delta, delta_ra))

        self.rotate(theta)

    def _align(self):
        """
        Aligns images via panning so their backgrounds stay consistent. Images requiring a pan greater than 1/2 the
         loaded image will be ignored.
        """
        x, y = WCS(self.header).all_world2pix(self.center[0], self.center[1], 0)

        if not(0 < x < self.get_data_size()[0] and 0 < y < self.get_data_size()[1]):
            self.logger.info("Pan out of range: ({}, {}) is greater than half the viewing window.".format(x, y))
        else:
            self.set_pan(x, y)

    def _key_press(self, canvas, keyname, opn, viewer):
        """
        Method called once a keyboard stoke has been detected. Using two un-bound keys, f & g, to cycle different
         cutout hdu's from the ObsRecord.
        Parameters canvas, opn, and viewer are all needed for the method to be called even though they are not
         directly used.

        :param canvas: Ginga DrawingCanvas Object
        :param keyname: Name of the key that has been pressed
        :param opn: str "key"
        :param viewer: Ginga EnhancedCanvasView object
        """
        self.logger.debug("Got key: {} from canvas: {} with opn: {} from viewer: {}".format(canvas,
                                                                                            keyname,
                                                                                            opn,
                                                                                            viewer))
        # Only step back if we aren't looking at a comparison images (as determined by the next_image keyword)
        if keyname == 'f':
            if self.next_image is not None:
                self.next_image = None
            else:
                self.obs_number -= 1
                key = self.downloader.image_key(self.candidate[self.obs_number])
                if key in self.comparison_images:
                    self.next_image = self.comparison_images[key]

        # Only step forward if this images doesn't have comparison image in the comparison image list.
        elif keyname == 'g':
            key = self.downloader.image_key(self.candidate[self.obs_number])
            if key in self.comparison_images and self.next_image is None:
                self.next_image = self.comparison_images[key]
            else:
                self.next_image = None
                self.obs_number += 1

        self.zoom = self.get_zoom()
        self.obs_number %= len(self.candidate)
        self._load()

    def clear_candidate_images(self):
        """
        Clear all the images associated with a candidate.
        """
        if self.candidate is None:
            return
        for obs_record in self.candidate:
            key = self.downloader.image_key(obs_record)
            if key in self.image_list:
                del(self.image_list[key])
            if key in self.astro_images:
                del(self.astro_images[key])
            if key in self.comparison_images:
                comp_key = self.comparison_images[key]
                if comp_key in self.comparison_images:
                    del(self.image_list[comp_key])
                if comp_key in self.astro_images:
                    del(self.astro_images[comp_key])

    def clear_viewer(self):
        """
        Clear the image in the viewer and any other objects drawn on the canvas.g
        """
        self.clear()
        self.canvas.delete_objects(self.canvas.get_objects())

    def buttons_on(self):
        """
        Activate most GUI buttons
        """
        self.next_set.set_enabled(True)
        self.previous_set.set_enabled(True)
        self.accept.set_enabled(True)
        self.reject.set_enabled(True)
        self.clear_button.set_enabled(True)
        self.load_json.set_enabled(True)
        self.reload_button.set_enabled(True)

    def buttons_off(self):
        """
        Deactivate some GUI buttons
        """
        self.next_set.set_enabled(False)
        self.previous_set.set_enabled(False)
        self.accept.set_enabled(False)
        self.reject.set_enabled(False)
        self.clear_button.set_enabled(False)
        self.load_json.set_enabled(False)
        self.reload_button.set_enabled(False)

    @property
    def center(self):
        """
        Returns the center of the image in ra/dec coordinates
        """
        if self._center is not None:
            return self._center

    @property
    def key(self):
        if self.next_image is not None:
            key = self.next_image
        else:
            key = self.downloader.image_key(self.candidate[self.obs_number])
        return key

    @property
    def loaded_hdu(self):
        """
        Return current HDU
        """
        # TODO: MEF
        key = self.key
        with self.lock:
            hdu = (isinstance(self.image_list[key], ApplyResult) and self.image_list[key].get()
                   or self.image_list[key])
            if isinstance(hdu, ApplyResult):
                self.logger.info("Loaded HDU is Apply result instance, not an HDU.")
                raise TypeError
            self.image_list[key] = hdu

        load_hdu = max_size = None
        for hdu in self.image_list[key]:
            if hdu.header['NAXIS'] == 0:
                continue
            size = hdu.header['NAXIS1'] * hdu.header['NAXIS2']
            if max_size is None or size > max_size:
                max_size = size
                load_hdu = hdu
        return load_hdu

    @property
    def header(self):
        """
        Return current HDU's header
        """
        return self.astro_images[self.key].get_header()

    @property
    def info(self):
        return "\n".join([x + " = " + str(self.header.get(x, "UNKNOWN")) for x in DISPLAY_KEYWORDS])


def main(params):

    ginga_logger = log.get_logger("ginga", options=params)

    if params.use_opencv:
        from ginga import trcalc
        try:
            trcalc.use('opencv')
        except Exception as ex:
            ginga_logger.warning("Error using OpenCL: {}".format(ex))

    if params.use_opencl:
        from ginga import trcalc
        try:
            trcalc.use('opencl')
        except Exception as ex:
            ginga_logger.warning("Error using OpenCL: {}".format(ex))

    app = Widgets.Application(logger=ginga_logger, host=params.host, port=params.port)

    #  create top level window
    window = app.make_window("Validate", wid='Validate')

    # our own viewer object, customized with methods (see above)
    ValidateGui(logging.getLogger('daomop'), window)

    try:
        app.start()

    except KeyboardInterrupt:
        ginga_logger.info("Terminating viewer...")
        window.close()
