from sensor_noise import *
import numpy as np
from math import *


def earth_sensor(true_value):
    """

    :param: true_value: true roll_angle of spacecraft, in radians
    :return: measured_value : roll_angle reading as measured by the earth_sensor, in radians

    """

    # values of Earth Sensor Characteristics taken from Sidi
    stat_noise_RMS_deg = radians(0.012)
    quantization_deg = radians(0.02)

    if true_value < 0:
        true_value_sign = -1
    else:
        true_value_sign = +1

    true_value_magnitude = abs(true_value)

    true_value_magnitude_deg = np.rad2deg(true_value_magnitude)

    multiple_number = np.floor(true_value_magnitude_deg/quantization_deg)
    remainder = true_value_magnitude_deg % quantization_deg
    if remainder >= 0.5*quantization_deg:
        multiple_number = multiple_number + 1
    else:
        pass

    measured_value_deg = true_value_sign * multiple_number * quantization_deg
    measured_value = np.deg2rad(measured_value_deg) + sensor_noise(stat_noise_RMS_deg)

    return measured_value
