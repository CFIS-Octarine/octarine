from astropy.time import Time

CFHT_QRUNS = {
    '17AQ01': (Time('2017-01-31 22:00:00'), Time('2017-02-03 22:00:00')),
    '17AQ03': (Time('2017-02-16 22:00:00'), Time('2017-03-04 22:00:00')),
    '17AQ06': (Time('2017-03-19 22:00:00'), Time('2017-04-04 22:00:00')),
    '17AQ08': (Time('2017-04-17 22:00:00'), Time('2017-05-04 22:00:00')),
    '17AQ10': (Time('2017-05-15 22:00:00'), Time('2017-05-31 22:00:00')),
    '17AQ13': (Time('2017-06-14 22:00:00'), Time('2017-06-28 22:00:00')),
    '17AQ16': (Time('2017-07-13 22:00:00'), Time('2017-08-01 22:00:00')),
    'default': (Time('2017-01-01 00:00:00'), Time('2021-01-01 00:00:00'))
}


def qrunid_start_date(qrunid):
    """

    :param qrunid: CFHT QRUN to give the start date of
    :return: start date of QRUN
    :rtype: Time
    """
    return CFHT_QRUNS.get(qrunid, CFHT_QRUNS['default'])[0].mjd


def qrunid_end_date(qrunid):
    """

    :param qrunid: CFHT QRUN to give the end date of
    :return: end date of QRUN
    :rtype: Time
    """
    return CFHT_QRUNS.get(qrunid, CFHT_QRUNS['default'])[1].mjd
