import numpy as np
import math


def sensor_noise(error_std_dev_deg):
    """
    :param: error_std_dev_deg : std. dev. of the statistical noise, in degrees
    :return: white_noise : statistical sensor noise, in radian
    """

    error_std_dev_rad = math.radians(error_std_dev_deg)
    white_noise = np.random.normal(0, error_std_dev_rad, 1)

    return white_noise
