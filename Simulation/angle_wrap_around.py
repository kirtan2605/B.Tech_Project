from math import *
import numpy as np


# Mika's Tech Blog
def truncated_remainder(dividend, divisor):
    divided_number = dividend / divisor
    divided_number = \
        -int(-divided_number) if divided_number < 0 else int(divided_number)
    truncatedremainder = dividend - divisor * divided_number
    return truncatedremainder


def transform_to_minus_pi_to_pi(input_angle):
    # revolutions = int((input_angle + np.sign(input_angle) * pi)/(2*pi))

    p1 = truncated_remainder(input_angle + np.sign(input_angle) * pi, 2 * pi)
    p2 = (np.sign(np.sign(input_angle) + 2 * (
                np.sign(fabs((truncated_remainder(input_angle + pi, 2 * pi)) / (2 * pi))) - 1))) * pi

    output_angle = p1 - p2

    return output_angle
