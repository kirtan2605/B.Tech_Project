import numpy as np


def PD_Control(ip, op_prev, Dt):

    if len(ip) < 2:
        ip = np.pad(ip, (2-len(ip), 0), 'constant', constant_values=(0, 0))

    # PD Control parameters
    kp = 0.00573
    kd = 0.3

    a = ((kd*2/Dt) + kp)
    b = (kp - (kd*2/Dt))

    # Bilinear Transform Implementation
    Controller_Output = -op_prev + a*ip[-1] + b*ip[-2]

    return Controller_Output
