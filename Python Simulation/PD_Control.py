import numpy as np


def PD_Control(ip, op_prev, Dt):

    if len(ip) < 3:
        ip = np.pad(ip, (3-len(ip), 0), 'constant', constant_values=(0, 0))

    #  PD Control parameters
    ki = 0
    kp = 0.00345
    kd = 10.6578

    # a = kp + kd + ki/2
    # b = -kp - 2*kd + ki/2
    # c = kd

    Controller_Output = kp*ip[-1] + kd*((ip[-1] - ip[-2])/Dt)

    return Controller_Output
