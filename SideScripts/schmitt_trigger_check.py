import matplotlib.pyplot as plt
from math import *
import numpy as np

def dual_schmitt_trigger(ip, op_prev):

    # Dual Schmitt Trigger Parameters
    OHP = 2.5
    OHM = -2.5
    OFF = 0

    ILM = -1.5
    IHM = -2.2
    ILP = 1.5
    IHP = 2.2

    if (ip < IHM ):
        pwpfm_op = OHM
    if (ip > IHP ):
        pwpfm_op = OHP
    if (ip > ILM and ip < ILP):
        pwpfm_op = OFF
    if (ip < ILM and ip > IHM) or (ip > ILP and ip < IHP) :
        pwpfm_op = op_prev

    return pwpfm_op

x = np.linspace(start = 0,stop = 20, num = 1001)

trigger_check_array = 4*np.sin(x) + 0.5*np.sin(50*x)

output = [0]

for i in range(1,len(trigger_check_array)) :
    output.append(dual_schmitt_trigger(trigger_check_array[i],output[-1]))


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
