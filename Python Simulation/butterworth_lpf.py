from scipy.signal import butter
from math import pi


def butterworth_lpf(f_s, f_c, order):
    # omega_c = 2 * pi * f_c  # Cut-off angular frequency
    omega_c = f_c              # Cut-off angular frequency
    omega_c_d = omega_c / f_s  # Normalized cut-off frequency (digital)

    # Design the digital Butterworth filter
    # the butter function expects the normalized frequency as a number from 0 to 1
    # so we divide it by pi before passing it to the butter function
    b, a = butter(order, omega_c_d / pi, btype='lowpass', analog=False, fs=f_s)

    return b, a
