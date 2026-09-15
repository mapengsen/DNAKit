"""Use the application's standard logging for imported scientific adapters."""

import logging


def get_logger(name):
    return logging.getLogger(name)
