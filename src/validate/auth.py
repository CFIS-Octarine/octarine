import getpass


def get_cadc_username():
    """
    get a target_name to use for locking and logging
    """
    return getpass.getuser()
