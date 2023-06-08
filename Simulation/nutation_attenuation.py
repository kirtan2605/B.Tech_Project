import numpy as np
from math import *
Dt = 0.01


thruster_on_counter = 0
thruster_off_counter = 0

inhibition = False

omega_nutation = 0.01414
    
#### add closed loop nutation frequency instead!!!!
print(omega_nutation)


f_nut = omega_nutation/(2*pi)
t_nut = 1/f_nut
inhibition_time = t_nut*0.5           # from Iwens
inhibition_step = int(round(inhibition_time / Dt))

####    Nutation Attentuation   ###
        # removing sudden sign change of control torque
        if (np.sign(control_torque_arr[i+1])*np.sign(control_torque_arr[i]) < 0) :
            control_torque_arr[i+1] = 0


        # nutation attenuation by seperating pules with half the nutation period.
        if not inhibition :
            if (control_torque_arr[i+1] != 0 and control_torque_arr[i] == 0) :
                thruster_on_step = i+1
                thruster_on_counter = thruster_on_counter + 1
                print("thruster on")
            if (control_torque_arr[i+1] == 0 and control_torque_arr[i] != 0) :
                thruster_off_step = i+1
                thruster_off_counter = thruster_off_counter + 1
                print("thruster off")

        if thruster_off_counter%2 != 0 :
            if i+1 < thruster_on_step + inhibition_step :
                inhibition = True
                control_torque_arr[i+1] = 0
            if i+1 == thruster_on_step + inhibition_step :
                print("inhibition ends")
                control_torque_arr[i+1] = control_torque_arr[i+1-inhibition_step]
                thruster_on_counter = thruster_on_counter + 1
            if i+1 > thruster_on_step + inhibition_step and i+1 < thruster_off_step + inhibition_step :
                control_torque_arr[i+1] = control_torque_arr[i+1-inhibition_step]
            if i+1 == thruster_off_step + inhibition_step :
                control_torque_arr[i+1] = control_torque_arr[i+1-inhibition_step]
                thruster_off_counter = thruster_off_counter + 1
                inhibition = False
