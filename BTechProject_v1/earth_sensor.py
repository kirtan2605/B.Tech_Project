from sensor_noise import *
import numpy as np


def earth_sensor(true_value):
    """

    :param: true_value: true roll_angle of spacecraft, in radians
    :return: measured_value : roll_angle reading as measured by the earth_sensor, in radians

    """

    # ....Sidi 8.3.2
    # modeling the earth sensor with order of accuracy : 0.02 deg = 0.000349 rad
    # this accuracy is accompanied by statistical noise : 0.03 deg RMS
    # modeling sensor noise as white noise with StdDev = RMS = 0.03 deg
    # ...DOUBT....

    stat_noise_RMS_deg = 0.03
    order_of_accuracy_deg = 0.02


    # check this snippet

    if true_value < 0:
        true_value_sign = -1
    else:
        true_value_sign = +1

    true_value_magnitude = abs(true_value)

    true_value_magnitude_deg = np.rad2deg(true_value_magnitude)

    multiple_number = np.floor(true_value_magnitude_deg/order_of_accuracy_deg)
    remainder = true_value_magnitude_deg % order_of_accuracy_deg
    if remainder >= 0.5*order_of_accuracy_deg:
        multiple_number = multiple_number + 1
    else:
        pass

    measured_value_deg = true_value_sign * multiple_number * order_of_accuracy_deg

    measured_value = np.deg2rad(measured_value_deg) + sensor_noise(stat_noise_RMS_deg)

    return measured_value
