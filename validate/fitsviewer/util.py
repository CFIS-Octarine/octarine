__author__ = "David Rusk <drusk@uvic.ca>"

from displayable import DisplayableImageSinglet


def display_hdulist(hdulist):
    DisplayableImageSinglet(hdulist).render()
