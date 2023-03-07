
def dual_schmitt_trigger(ip, op_prev):

    # Dual Schmitt Trigger Parameters calculation!!
    
    Um = 0.05
    OFF = 0

    h = 0.045
    U_on = 0.05
    U_off = (U_on - h)

    # returning pwpf_output
    if (ip < -U_on ):
        pwpfm_op = -Um
    if (ip > U_on ):
        pwpfm_op = Um
    if (ip > -U_off and ip < U_off):
        pwpfm_op = OFF
    if (ip < -U_off and ip > -U_on) or (ip > U_off and ip < U_on) :
        pwpfm_op = op_prev

    return pwpfm_op