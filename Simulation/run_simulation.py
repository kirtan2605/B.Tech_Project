import matplotlib.pyplot as plt

from angle_wrap_around import *
from PD_Control import *
from disturbance_torque import *
from earth_sensor import *

from butterworth_lpf import *
from low_pass_filter import *
from schmitt_trigger import *
from pwpfm import *

def plot_yaw_vs_roll(x_arr):
    # plotting the result
    plt.plot(np.rad2deg(x_arr[2, :]), np.rad2deg(x_arr[0, :]),linewidth=1)
    # plt.plot..... plot another data in same plot if needed
    plt.title('Roll Angle vs Yaw Angle', fontsize=12)
    plt.ylabel('Roll Angle (degrees)', fontsize=12)
    plt.xlabel('Yaw Angle (degrees)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    # plt.axis
    # plt.legend()
    plt.show()



def run_simulation(roll_desired, alpha_d, runtime_parameters, system_variables, sampling_parameters, initial_conditions):
    t_start, t_end, Dt = runtime_parameters
    Ix, Iy, Iz, omega0, OmegaN, h = system_variables
    output_time, one = sampling_parameters
    T_sam = output_time

    a = 4 * (omega0 ** 2) * (Iy - Iz)
    b = -1 * omega0 * (Ix - Iy + Iz)
    c = (omega0 ** 2) * (Iy - Ix)

    # x_dot = A.x + B.Md + B.Mc
    B = [[0, 0], [1/Ix, 0], [0, 0], [0, 1/Iz]]
    A1 = (-1 * (a + omega0 * h)) / Ix
    A2 = (-(b + h)) / Ix
    A3 = (b + h) / Iz
    A4 = (-1 * (c + omega0 * h)) / Iz
    A = [[0, 1, 0, 0], [A1, 0, 0, A2], [0, 0, 0, 1], [0, A3, A4, 0]]

    # calculating number of timestamps in simulation
    n_steps = int(round(t_end - t_start) / Dt)
    t_arr = np.append((np.arange(t_start, t_end, Dt)), t_end)

    # x = [phi, phi_dot, psi, psi_dot]
    x_arr = np.zeros((4, n_steps + 1))

    # external torque vector M_arr = [Mdx, Mdz, Mcx, Mcz]
    Md_arr = np.array(disturbance_torques(n_steps + 1))
    Mc_arr = np.array(np.zeros((2, n_steps + 1)))

    # The Earth Sensor measures the pitch and roll angles only
    # the roll rate is  calculated from consecutive roll measurements
    control_torque_arr = np.zeros(n_steps + 1)
    positive_control_torque_arr = np.zeros(n_steps + 1)
    negative_control_torque_arr = np.zeros(n_steps + 1)
    roll_error_arr = np.zeros(n_steps + 1)
    phi_measured_arr = np.zeros(n_steps + 1)
    phi_measured_lpf_arr = np.zeros(n_steps + 1)
    roll_error_measured_arr = np.zeros(n_steps + 1)
    
    pwpfm_error_arr = np.zeros(n_steps + 1)
    pwpfm_error_lpf_arr = np.zeros(n_steps + 1)

    # adding initial conditions
    x_arr[:, 0] = initial_conditions
    roll_error_arr[0] = roll_desired - x_arr[0, 0]

    # Earth Sensor and Controller Delay Implementation
    output_step = int(output_time / Dt)


    T_c = 1  # Torque Magnitude Multiplier

    Wn_nut = 0.006688126294945937           # nutation frequency in rad/sec
    f_nut = Wn_nut/(2*pi)
    t_nut = 1/f_nut
    inhibition_time = t_nut/2
    inhibition_step = int(round(inhibition_time / Dt))
    print("inhibition steps : ", inhibition_step)
    print("inhibition time : ", inhibition_time)

    # using a Butterworth Filter of order 1 to filter out the sensor noise
    filter_order = 1  # Order of the butterworth filter
    f_sample = 1 / T_sam  # Sample frequency in Hz
    f_cutoff_rps = Wn_nut * 10  # Cut-off frequency in rad/sec (closed loop wn*7.5)
    b, a = butterworth_lpf(f_sample, f_cutoff_rps, filter_order)

    Tc_controller_output = 0

    #positive_thruster_on_counter = 0
    #positive_thruster_off_counter = 0

    thruster_on_counter = 0
    thruster_off_counter = 0

    #negative_thruster_on_counter = 0
    #negative_thruster_off_counter = 0

    #positive_inhibition = False
    #negative_inhibition = False

    inhibition = False

    for i in range(0, n_steps):

        x = x_arr[:, i]
        roll_error_arr[i] = roll_desired - x_arr[0, i]

        # modelling of sample-and-hold earth sensor
        if i % output_step == 0:
            phi_measured_arr[i] = earth_sensor(x_arr[0, i])
        else:
            phi_measured_arr[i] = phi_measured_arr[i - 1]

        # applying low pass filter to filter out noise
        if i - filter_order < 0:
            phi_measured_lpf_arr[i] = low_pass_filter(phi_measured_lpf_arr[0:i + 1], phi_measured_arr[0:i + 1], b, a,
                                                      filter_order)
        else:
            phi_measured_lpf_arr[i] = low_pass_filter(phi_measured_lpf_arr[i - filter_order:i + 1],
                                                      phi_measured_arr[i - filter_order:i + 1], b, a, filter_order)

        roll_error_measured_arr[i] = roll_desired - phi_measured_lpf_arr[i]

        # implementing roll dead-band # check
        # if abs(roll_error_measured_arr[i]) < roll_deadband_rad:
        #    roll_error_measured_arr[i] = 0

        # calculating rate of change of state variables
        x_dot = np.matmul(A, x) + np.matmul(B, Mc_arr[:, i]) + np.matmul(B, Md_arr[:, i])

        # updating state variables
        x_arr[:, i + 1] = x + np.dot(x_dot, Dt)
        x_arr[0, i + 1] = transform_to_minus_pi_to_pi(x_arr[0, i + 1])
        x_arr[2, i + 1] = transform_to_minus_pi_to_pi(x_arr[2, i + 1])

        if i < 2:
            Tc_controller_output = PD_Control(roll_error_measured_arr[0:i], Tc_controller_output, Dt)
        else:
            Tc_controller_output = PD_Control(roll_error_measured_arr[i-2:i + 1], Tc_controller_output, Dt)

        control_torque_arr[i+1], pwpfm_error_arr[i+1], pwpfm_error_lpf_arr[i+1] =  pwpfm(T_c*Tc_controller_output, pwpfm_error_arr[i], pwpfm_error_lpf_arr[i], control_torque_arr[i], T_sam)

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




        '''

        # nutation attenuation by seperating pules with half the nutation period.
        if not positive_inhibition :
            if (control_torque_arr[i+1] > 0 and control_torque_arr[i] <= 0) :
                positive_thruster_on_step = i+1
                positive_thruster_on_counter = positive_thruster_on_counter + 1
                print("positive thruster on")
            if (control_torque_arr[i+1] <= 0 and control_torque_arr[i] > 0) :
                positive_thruster_off_step = i+1
                positive_thruster_off_counter = positive_thruster_off_counter + 1
                print("positive thruster off")

        if positive_thruster_off_counter%2 != 0 :
            if i+1 < positive_thruster_on_step + inhibition_step :
                positive_inhibition = True
                control_torque_arr[i+1] = 0
            if i+1 == positive_thruster_on_step + inhibition_step :
                print("positive inhibition ends")
                control_torque_arr[i+1] = control_torque_arr[i+1-inhibition_step]
                positive_thruster_on_counter = positive_thruster_on_counter + 1
            if i+1 > positive_thruster_on_step + inhibition_step and i+1 < positive_thruster_off_step + inhibition_step :
                control_torque_arr[i+1] = control_torque_arr[i+1-inhibition_step]
            if i+1 == positive_thruster_off_step + inhibition_step :
                control_torque_arr[i+1] = control_torque_arr[i+1-inhibition_step]
                positive_thruster_off_counter = positive_thruster_off_counter + 1
                positive_inhibition = False
           
        if not negative_inhibition :
            if (control_torque_arr[i+1] > 0 and control_torque_arr[i] <= 0) :
                negative_thruster_on_step = i+1
                negative_thruster_on_counter = negative_thruster_on_counter + 1
                print("negative thruster on")
            if (control_torque_arr[i+1] <= 0 and control_torque_arr[i] > 0) :
                negative_thruster_off_step = i+1
                negative_thruster_off_counter = negative_thruster_off_counter + 1
                print("negative thruster off")

        if negative_thruster_off_counter%2 != 0 :
            if i+1 < negative_thruster_on_step + inhibition_step :
                negative_inhibition = True
                control_torque_arr[i+1] = 0
            if i+1 == negative_thruster_on_step + inhibition_step :
                print("negative inhibition ends")
                control_torque_arr[i+1] = control_torque_arr[i+1-inhibition_step]
                negative_thruster_on_counter = negative_thruster_on_counter + 1
            if i+1 > negative_thruster_on_step + inhibition_step and i+1 < negative_thruster_off_step + inhibition_step :
                control_torque_arr[i+1] = control_torque_arr[i+1-inhibition_step]
            if i+1 == negative_thruster_off_step + inhibition_step :
                control_torque_arr[i+1] = control_torque_arr[i+1-inhibition_step]
                negative_thruster_off_counter = negative_thruster_off_counter + 1
                negative_inhibition = False
        '''

        # effects of actuation set-up i.e. offset nature of thrusters
        offset_actuation = [cos(radians(alpha_d)), -1 * sin(radians(alpha_d))]
        Mc_arr[:, i + 1] = np.multiply(offset_actuation, control_torque_arr[i + 1])


    # plotting the result
    plt.plot(t_arr, np.rad2deg(x_arr[0, :]), linewidth=1, label='Roll Angle (degrees)')
    # plt.plot(t_arr, np.rad2deg(phi_measured_arr), linewidth=0.5, label='Measured Roll')
    plt.plot(t_arr, np.rad2deg(phi_measured_lpf_arr), linewidth=1, label='Measured Roll LPF')
    # plt.plot..... plot another data in same plot if needed
    plt.title('Roll Angle vs time', fontsize=12)
    plt.xlabel('t (seconds)', fontsize=12)
    plt.ylabel('Angle (degrees)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    # plt.axis
    plt.legend()
    plt.show()

    # plotting the result
    plt.plot(t_arr, np.rad2deg(x_arr[2, :]), linewidth=1, label='Yaw Angle (degrees)')
    # plt.plot..... plot another data in same plot if needed
    plt.title('Yaw Angle vs time', fontsize=12)
    plt.xlabel('t (seconds)', fontsize=12)
    plt.ylabel('Angle (degrees)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    # plt.axis
    # plt.legend()
    plt.show()

    # plotting the result
    plt.plot(t_arr, control_torque_arr[:],linewidth=1, label='Control Torque')
    #plt.plot(t_arr, control_torque_arr[:],marker = '.',linewidth=1, label='Control Torque')
    # plt.plot..... plot another data in same plot if needed
    plt.title('Control Torque vs time', fontsize=12)
    plt.xlabel('t (seconds)', fontsize=12)
    plt.ylabel('Torque (Nm)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    # plt.axis
    # plt.legend()
    plt.show()

    # plot_yaw_vs_roll(x_arr)
