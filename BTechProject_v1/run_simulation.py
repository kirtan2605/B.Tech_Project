
import matplotlib.pyplot as plt

from angle_wrap_around import *
from PDF_Control import *
from disturbance_torque import *
from earth_sensor import *



def run_simulation(roll_desired, alpha_d, runtime_parameters, system_variables, initial_conditions):

    t_start, t_end, Dt = runtime_parameters
    Ix, Iy, Iz, omega0, h = system_variables

    a = 4 * (omega0 ** 2) * (Iy - Iz)
    b = -1 * omega0 * (Ix - Iy + Iz)
    c = (omega0 ** 2) * (Iy - Ix)

    # calculating number of timestamps in simulation
    n_steps = int(round(t_end - t_start) / Dt)

    t_arr = np.append((np.arange(t_start, t_end, Dt)), t_end)

    # x = [phi, phi_dot, psi, psi_dot]
    x_arr = np.zeros((4, n_steps + 1))

    # external torque vector M_arr = [Mdx, Mdz, Mcx, Mcz]
    M_arr = np.concatenate((disturbance_torques(n_steps+1), np.zeros((2, n_steps+1))))

    # x_dot = A.x + B.U

    B = [[0, 0, 0, 0], [1 / Ix, 0, 1 / Ix, 0], [0, 0, 0, 0], [0, 1 / Iz, 0, 1 / Iz]]

    A1 = (-1 * (a + omega0 * h)) / Ix
    A2 = (-(b + h)) / Ix
    A3 = (b + h) / Iz
    A4 = (-1 * (c + omega0 * h)) / Iz
    A = [[0, 1, 0, 0], [A1, 0, 0, A2], [0, 0, 0, 1], [0, A3, A4, 0]]


    # The Earth Sensor measures the pitch and roll angles only
    # the phi_dot_measured is actually calculated from consecutive roll measurements

    roll_error_arr = np.zeros(n_steps + 1)
    phi_measured_arr = np.zeros(n_steps + 1)
    roll_error_measured_arr = np.zeros(n_steps + 1)



    # adding initial conditions
    x_arr[:, 0] = initial_conditions
    roll_error_arr[0] = roll_desired - x_arr[0, 0]
    phi_measured_arr[0] = earth_sensor(x_arr[0, 0])
    roll_error_measured_arr[0] = roll_desired - phi_measured_arr[0]


    # let's say that the measurement is taken by the Earth Sensor every t_measure seconds
    # t_measure must be greater than simulation Dt, hence let
    # t_measure = n * Dt        n > 1
    # output_rate = 1/t_measure
    # output_step = n
    # sample taken at every nth step of the simulation iterate
    # we assume sampling always happens from the first step of iteration
    # try to keep sampling time such that it is an integral multiple of Dt
    output_time = 0.2               # time between each measurement of earth sensor
    output_rate = 1/output_time
    output_step = int(output_time/Dt)

    print("\nEarth Sensor Output Rate : ", output_rate, "Hz")
    print("No Earth Sensor Information is lost due to sampling")
    print("since Earth Sensor Output Rate <= Sampling Frequency\n")


    for i in range(1, n_steps + 1):
        x = x_arr[:, i-1]

        if i == 1:
            DrollerrorDt = 0
        else:
            DrollerrorDt = (roll_error_measured_arr[i-1] - roll_error_measured_arr[i-2]) / Dt

        # calculating control torque
        Tc_controller_output = PDF_Control_Dougherty(roll_error_measured_arr[i-1], DrollerrorDt)

        Tc = Tc_controller_output

        # effects of actuation set-up i.e. offset nature of thrusters
        offset_actuation = [cos(radians(alpha_d)), -1 * sin(radians(alpha_d))]
        M_arr[2:, i] = np.multiply(offset_actuation, Tc)

        # external torques
        U = M_arr[:, i]

        # calculating rate of change of state variables
        x_dot = np.matmul(A, x) + np.matmul(B, U)

        # updating state variables
        x_arr[:, i] = x + np.dot(x_dot, Dt)
        x_arr[0, i] = transform_to_minus_pi_to_pi(x_arr[0, i])
        x_arr[2, i] = transform_to_minus_pi_to_pi(x_arr[2, i])

        roll_error_arr[i] = roll_desired - x_arr[0, i]

        if (i-1) % output_step == 0:
            phi_measured_arr[i] = earth_sensor(x_arr[0, i])
        else:
            phi_measured_arr[i] = phi_measured_arr[i-1]

        roll_error_measured_arr[i] = roll_desired - phi_measured_arr[i]

    # plotting the result
    plt.plot(t_arr, np.rad2deg(x_arr[0, :]), linewidth=1, label='Roll Angle (degrees)')
    plt.plot(t_arr, np.rad2deg(phi_measured_arr), linewidth=1, label='Measured Roll')
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
    plt.legend()
    plt.show()