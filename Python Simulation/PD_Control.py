import numpy as np


def PD_Control(ip, op_prev, Dt):

    if len(ip) < 2:
        ip = np.pad(ip, (2-len(ip), 0), 'constant', constant_values=(0, 0))

    #  PD Control parameters
    ki = 0
    kp = 0.00345
    kd = 10.6578

    ## simple implementation
    # Controller_Output = kp*ip[-1] + kd*((ip[-1] - ip[-2])/Dt)

    a = ((kd*2/Dt) + kp)
    b = (kp - (kd*2/Dt))

    ## Bilinear Transform Implementation
    Controller_Output = -op_prev + a*ip[-1] + b*ip[-2]

    return Controller_Output
