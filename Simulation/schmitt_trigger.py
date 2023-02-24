
def dual_schmitt_trigger(ip, op_prev):

    # Dual Schmitt Trigger Parameters
    OHP = 0.01
    OHM = -0.01
    OFF = 0

    ILM = -0.001
    IHM = -0.03
    ILP = 0.001
    IHP = 0.03

    # returning pwpf_output
    if (ip < IHM ):
        pwpfm_op = OHM
    if (ip > IHP ):
        pwpfm_op = OHP
    if (ip > ILM and ip < ILP):
        pwpfm_op = OFF
    if (ip < ILM and ip > IHM) or (ip > ILP and ip < IHP) :
        pwpfm_op = op_prev

    return pwpfm_op