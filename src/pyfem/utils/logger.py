import logging


def set_logger(props):
    level = "normal"

    if hasattr(props, "logger"):
        level = props.logger.level

        if level not in ["normal", "info", "debug", "critical", "warning", "silent"]:
            raise NotImplementedError(
                'Logger level should be "normal", "info", "debug", "critical", "silent" or "warning"')

    logger = logging.getLogger()
    handler = logging.StreamHandler()

    if level == "debug":
        formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s')
        logger.setLevel(logging.DEBUG)
    elif level == "critical" or level == "silent":
        formatter = logging.Formatter('  %(message)s')
        logger.setLevel(logging.CRITICAL)
    elif level == "warning":
        formatter = logging.Formatter('  %(message)s')
        logger.setLevel(logging.WARNING)
    else:
        formatter = logging.Formatter('  %(message)s')
        logger.setLevel(logging.INFO)

    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger


def get_logger():
    return logging.getLogger()
