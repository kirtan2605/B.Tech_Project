import numpy as np


def PD_Control(ip, op_prev, Dt):

    if len(ip) < 2:
        ip = np.pad(ip, (2-len(ip), 0), 'constant', constant_values=(0, 0))

    # PD Control parameters
    kp = 0.004275
    kd = 0.673


    a = ((kd*2/Dt) + kp)
    b = (kp - (kd*2/Dt))

    # Bilinear Transform Implementation
    Controller_Output = + a*ip[-1] + b*ip[-2] - op_prev 

    return Controller_Output
