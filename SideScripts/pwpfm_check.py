import matplotlib.pyplot as plt
from math import *
import numpy as np


def dual_schmitt_trigger(ip, op_prev):

    # Dual Schmitt Trigger Parameters calculation!!
    
    Um = 2.5
    OFF = 0

    h = 0.7
    U_on = 2.2
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


n = 1000
t_start = 0
t_stop = 10
Dt = (t_stop - t_start)/n

pwpfm_error_arr = np.zeros(n+2)
pwpfm_error_lpf_arr = np.zeros(n+2)

x = np.linspace(start = t_start,stop = t_stop, num = n+1)

trigger_check_array = 4*np.sin(x) + 0.5*np.sin(50*x)
A = 0.5
B = 0.3
C = 0.25
trigger_check_array = 0.25*np.sin(A*x**2 + B*x + C) + 0*np.cos(50*x)

output = np.zeros(n+2)

for i in range(0,len(trigger_check_array)) :
    output[i+1], pwpfm_error_arr[i+1], pwpfm_error_lpf_arr[i+1]  = pwpfm(trigger_check_array[i], pwpfm_error_arr[i], pwpfm_error_lpf_arr[i], output[i], Dt)


Um = 2.5
OFF = 0

h = 0.7
U_on = 2.2
U_off = (U_on - h)

# plotting the result
plt.plot(output, linewidth=1, label='Output')
plt.plot(trigger_check_array, linewidth=1, label='Input')
plt.title('Schmitt Trigger Output Check', fontsize=12)
plt.xlabel('Input Iteration', fontsize=12)
plt.ylabel('Input / Output', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.legend()
plt.show()
