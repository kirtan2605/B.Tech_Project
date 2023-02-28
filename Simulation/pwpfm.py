from schmitt_trigger import *

def pwpfm(ip, error_prev, error_lpf_prev,  op_prev, Dt):

    # pwpfm gain parameter
    k_pm = 20

    # lpf parameters
    km = 1
    tm = 0.1
    
    ip = k_pm*ip

    error = ip - op_prev

    T = (2*tm)/Dt
    error_lpf = (km*(error + error_prev) + T*error_lpf_prev)/(T+1)

    pwpfm_value = dual_schmitt_trigger(error_lpf, op_prev)

    #pwm_value = 0

    return_list = [pwpfm_value, error, error_lpf]
    return return_list
