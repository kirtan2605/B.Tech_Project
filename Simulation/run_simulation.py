import matplotlib.pyplot as plt

from angle_wrap_around import *
from PD_Control import *
from disturbance_torque import *
from earth_sensor import *

from butterworth_lpf import *
from low_pass_filter import *
from schmitt_trigger import *
from pwpfm import *
from solver import *

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



def run_simulation(roll_desired, alpha_d, runtime_parameters, sampling_parameters, initial_conditions):
    t_start, t_end, Dt = runtime_parameters
    sampling_time, controller_time = sampling_parameters
    
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
    sampling_step = int(sampling_time / Dt)
    #controller_step = int(controller_time / Dt)


    T_c = 100  # Torque Magnitude Multiplier

    omega_nutation = 0.01414        # closed loop nutation frequency

    # using a Butterworth Filter of order 1 to filter out the sensor noise
    filter_order = 1  # Order of the butterworth filter
    f_sample = 1 / sampling_time  # Sample frequency in Hz
    f_cutoff_rps = omega_nutation * 10  # Cut-off frequency in rad/sec (closed loop wn*7.5)
    b, a = butterworth_lpf(f_sample, f_cutoff_rps, filter_order)

    Tc_controller_output = 0


    for i in range(0, n_steps):

        x = x_arr[:, i]
        roll_error_arr[i] = roll_desired - x_arr[0, i]

        # modelling of sample-and-hold earth sensor
        if i % sampling_step == 0:
            phi_measured_arr[i] = earth_sensor(x_arr[0, i])
        else:
            phi_measured_arr[i] = phi_measured_arr[i-1]

        # applying low pass filter to filter out noise
        if i - filter_order < 0:
            phi_measured_lpf_arr[i] = low_pass_filter(phi_measured_lpf_arr[0:i+1], phi_measured_arr[0:i+1], b, a,
                                                      filter_order)
        else:
            phi_measured_lpf_arr[i] = low_pass_filter(phi_measured_lpf_arr[i-filter_order:i+1],
                                                      phi_measured_arr[i-filter_order:i+1], b, a, filter_order)

        #phi_measured_lpf_arr[i] = phi_measured_arr[i]

        roll_error_measured_arr[i] = roll_desired - phi_measured_lpf_arr[i]

        x_arr[:,i+1] = rk4(ode_system, i, x, Mc_arr[:, i], Md_arr[:, i], Dt)

        x_arr[0, i+1] = transform_to_minus_pi_to_pi(x_arr[0, i+1])
        x_arr[2, i+1] = transform_to_minus_pi_to_pi(x_arr[2, i+1])


        Tc_controller_output = PD_Control(roll_error_measured_arr[i-2:i+1], Tc_controller_output, sampling_time)
        control_torque_arr[i+1] = T_c*Tc_controller_output
        
        
        control_torque_arr[i+1], pwpfm_error_arr[i+1], pwpfm_error_lpf_arr[i+1] =  pwpfm(T_c*Tc_controller_output, pwpfm_error_arr[i], pwpfm_error_lpf_arr[i], control_torque_arr[i], sampling_time)
    
        
        # effects of actuation set-up i.e. offset nature of thrusters
        offset_actuation = [cos(radians(alpha_d)), -1 * sin(radians(alpha_d))]
        Mc_arr[:, i + 1] = np.multiply(offset_actuation, control_torque_arr[i + 1])


    
    # plotting the result
    plt.plot(t_arr/60, np.rad2deg(x_arr[0, :]), linewidth=1, label='Roll Angle (degrees)')
    #plt.plot(t_arr/60, np.rad2deg(phi_measured_arr), linewidth=0.5, label='Measured Roll')
    #plt.plot(t_arr/60, np.rad2deg(phi_measured_lpf_arr), linewidth=1, label='Measured Roll LPF')
    plt.axhline(y=degrees(roll_desired), color='r', linestyle='-')
    plt.axhline(y=degrees(roll_desired) + 0.05, color='r', linestyle=':')
    plt.axhline(y=degrees(roll_desired) - 0.05, color='r', linestyle=':')
    # plt.plot..... plot another data in same plot if needed
    plt.title('Roll Angle vs time', fontsize=12)
    plt.xlabel('t (minutes)', fontsize=12)
    plt.ylabel('Angle (degrees)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    # plt.axis
    plt.legend()
    plt.show()


    # plotting the result
    plt.plot(np.rad2deg(x_arr[0, :]),np.rad2deg(x_arr[1, :]), linewidth=1, label='Roll')
    plt.plot(np.rad2deg(x_arr[2, :]),np.rad2deg(x_arr[3, :]), linewidth=1, label='Yaw')
    #plt.axhline(y=degrees(roll_desired), color='r', linestyle='-')
    #plt.axhline(y=degrees(roll_desired) + 0.05, color='r', linestyle=':')
    #plt.axhline(y=degrees(roll_desired) - 0.05, color='r', linestyle=':')
    # plt.plot..... plot another data in same plot if needed
    plt.title('Rate vs Angle', fontsize=12)
    plt.xlabel('Angle (deg)', fontsize=12)
    plt.ylabel('Rate (deg/s)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    # plt.axis
    plt.legend()
    plt.show()

    # plotting the result
    plt.plot(t_arr/60, np.rad2deg(x_arr[2, :]), linewidth=1, label='Yaw Angle (degrees)')
    # plt.plot..... plot another data in same plot if needed
    plt.axhline(y=0, color='r', linestyle='-')
    plt.axhline(y= 0.35, color='r', linestyle=':')
    plt.axhline(y= -0.35, color='r', linestyle=':')
    plt.title('Yaw Angle vs time', fontsize=12)
    plt.xlabel('t (minutes)', fontsize=12)
    plt.ylabel('Angle (degrees)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    # plt.axis
    # plt.legend()
    plt.show()
    

    # plotting the result
    plt.plot(t_arr/60, np.rad2deg(roll_error_measured_arr[:]), linewidth=1, label='Roll Error Measured (degrees)')
    # plt.plot..... plot another data in same plot if needed
    plt.title('Roll Error Measured vs time', fontsize=12)
    plt.xlabel('t (minutes)', fontsize=12)
    plt.ylabel('Angle (degrees)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    # plt.axis
    plt.legend()
    plt.show()


    # plotting the result
    plt.plot(t_arr/60, control_torque_arr[:],linewidth=1, label='Control Torque')
    #plt.plot(t_arr, control_torque_arr[:],marker = '.',linewidth=1, label='Control Torque')
    # plt.plot..... plot another data in same plot if needed
    plt.title('Control Torque vs time', fontsize=12)
    plt.xlabel('t (minutes)', fontsize=12)
    plt.ylabel('Torque (Nm)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    # plt.axis
    # plt.legend()
    plt.show()

    # plot_yaw_vs_roll(x_arr)

    # plotting the result
    plt.plot(np.rad2deg(x_arr[2, :]), np.rad2deg(x_arr[0, :]), linewidth=1, label='Roll Angle (degrees)')
    # plt.plot..... plot another data in same plot if needed
    plt.title('Roll Angle vs Yaw Angle', fontsize=12)
    plt.xlabel(' Yaw Angle (degreed))', fontsize=12)
    plt.ylabel(' Roll Angle (degrees)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    # plt.axis
    plt.legend()
    plt.show()