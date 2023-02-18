from scipy.signal import butter
from math import pi


def butterworth_lpf(f_s, f_c, order):
    # input cut-off frequency is already in angular frequency
    omega_c = f_c              # Cut-off angular frequency
    omega_c_d = omega_c / f_s  # Normalized cut-off frequency (digital)

    b, a = butter(order, omega_c_d/pi, btype='lowpass', analog=False, fs=f_s)

    return b, a
