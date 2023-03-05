from math import *

def system_variables():
    # calculating spacecraft moments of inertia
    # (INSAT-3DR : [2211, 2.4, 1.6, 1.5])
    mass = 2211             # mass in kg
    height = 2.4            # height in m
    width = 1.6             # width in m
    depth = 1.5             # depth in m
    # by convention height > width > depth
    Iz = mass*(height**2 + width**2)/12
    Ix = mass*(height**2 + depth**2)/12
    Iy = mass*(width**2 + depth**2)/12

    ## defining the moments of inertia
    # Ix = 800
    # Iy = 680
    # Iz = 1000

    moments_of_inertia = (Ix, Iy, Iz)

    return moments_of_inertia

def calculate_parameters(a_assumed):

    (Ix, Iy, Iz) = system_variables()

    wo = 2*pi/86400         # orbit frequency in rad/sec
    wo = 7.236e-05
    xi = 0.7                # damping coefficient of closed loop poles

    Tdx_max = 5e-6          # maximum magnitude of disturbance torque
    Tdz_max = 5e-6          # maximum magnitude of disturbance torque
    phi_ss = 0.05          # steady state error in roll in Deg
    psi_ss = 0.4           # steady state error in yaw in Deg

    phi_ss = radians(phi_ss)
    psi_ss = radians(psi_ss)


    # calculating h, kx without the approximation kx >> wo*h !!
    Kx = (Tdx_max*(psi_ss/phi_ss) - Tdz_max)/(psi_ss - a_assumed*phi_ss)
    h = (Tdx_max/phi_ss - Kx)/wo

    # calculating kxd,a_calculated, wn1, wn2
    A = sqrt(((wo * wo * h * h) + (wo * h * Kx)) / (Ix * Iz))
    B = 1 / (2 * Ix * xi)
    C = (wo * h * (Ix + Iz) + Iz * Kx + h * h) / (Ix * Iz)

    Kxd = sqrt(((C + A * (2 - 4 * xi * xi)) * Kx) / ((Kx * B * B) - (2 * xi * A * B) + (h * wo) / (Ix * Iz)))
    a_psi = (2 * xi * A * B * Kxd * Ix * Iz - wo * Kxd * h) / (h * Kx)
    Wn1 = ((Kxd * B) + sqrt(Kxd * Kxd * B * B - 4 * A)) / 2
    Wn2 = A / Wn1

    calculated_parameters = (wo, xi, h, Kx, Kxd, a_psi, Wn1, Wn2)

    return calculated_parameters

def display_parameters(parameters):
    sys_var = system_variables()
    print("\nSystem Parameters")
    print("Ix :", sys_var[0])
    print("Iy :", sys_var[1])
    print("Iz :", sys_var[2] , "\n")

    print("\nSimulation Parameters")
    print("w0 :", parameters[0])
    print("xi :", parameters[1])
    print("h :", parameters[2])
    print("Kx: ", parameters[3])
    print("Kxd: ", parameters[4])
    print("a_psi: ", parameters[5])
    print("alpha: ", degrees(atan(parameters[5])))
    print("Wn1: ", parameters[6])
    print("Wn2: ", parameters[7] , "\n")

# set error tolerence in alpha
alpha_tolerance_deg = 0.001            # tolerance in alpha in Deg
alpha_tolerance = radians(alpha_tolerance_deg)

alpha_guess_deg = 5
alpha_guess = radians(alpha_guess_deg)
a_guess = tan(alpha_guess)

simulation_parameters = calculate_parameters(a_guess)
a_calculated = simulation_parameters[3]
alpha_calculated = atan(a_calculated)

alpha_error = abs(alpha_calculated - alpha_guess)

while alpha_error > alpha_tolerance:

    alpha_guess = alpha_guess - (alpha_guess - alpha_calculated)/2
    a_guess = tan(alpha_guess)

    simulation_parameters = calculate_parameters(a_guess)
    a_calculated = simulation_parameters[3]
    alpha_calculated = atan(a_calculated)

    alpha_error = abs(alpha_calculated - alpha_guess)

display_parameters(simulation_parameters)