import logging
import os

def init_logger(logger_name="pfp_log", log_file_name="pfp.log", to_file=True, to_console=True):
    """
    Purpose:
     Returns a logger object.
    Usage:
     logger = pfp_log.init_logger()
    Author: PRI with acknowledgement to James Cleverly
    Date: September 2016
    """
    # create the logfiles directory if it does not exist
    if not os.path.exists("logfiles"):
        os.makedirs("logfiles")
    # create formatter
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s','%H:%M:%S')
    # create the logger
    logger = logging.getLogger(name=logger_name)
    logger.setLevel(logging.DEBUG)
    # create file handler if requested
    if to_file:
        log_file_path = os.path.join("logfiles", log_file_name)
        fh1 = logging.FileHandler(log_file_path)
        fh1.setLevel(logging.DEBUG)
        # add formatter to file handler
        fh1.setFormatter(formatter)
        # add the file handler to the logger
        logger.addHandler(fh1)
        # set up a separate file for errors
        error_file_name = log_file_name.replace(".log","_errors.log")
        error_file_path = os.path.join("logfiles", error_file_name)
        fh2 = logging.FileHandler(error_file_path)
        fh2.setLevel(logging.ERROR)
        fh2.setFormatter(formatter)
        logger.addHandler(fh2)
    # create console handler if requested
    if to_console:
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        # add formatter to console handler
        ch.setFormatter(formatter)
        # add the console handler to the logger
        logger.addHandler(ch)

    return logger
