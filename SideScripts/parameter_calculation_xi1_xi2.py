from math import *
import numpy as np

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

def display_parameters(parameters):
    Ix, Iy, Iz = system_variables()
    print("\nSystem Parameters")
    print("Ix :", Ix)
    print("Iy :", Iy)
    print("Iz :", Iz , "\n")

    wo, xi1, xi2, h, Kx, Kxd, a_psi, Wn1, Wn2 = parameters
    print("\nSimulation Parameters")
    print("w0 :", wo)
    print("xi1 :", xi1)
    print("xi2 :", xi2)
    print("h :", h)
    print("Kx: ", Kx)
    print("Kxd: ", Kxd)
    print("a_psi: ", a_psi)
    print("alpha: ", degrees(atan(a_psi)))
    print("Wn1: ", Wn1)
    print("Wn2: ", Wn2, "\n")

     # verify if the equations hold for epsilon e
    e = 0.01
    RHS1 = (Iz*Kxd)
    RHS2 = (wo*h*(Ix+Iz) + Iz*Kx + h*h + a_psi*h*Kxd)
    RHS3 = (a_psi*h*Kx + wo*h*Kxd)
    RHS4 = (wo*wo*h*h + wo*h*Kx)
    LHS1 = (2*(xi1*Wn1 + xi2*Wn2))*(Ix*Iz)
    LHS2 = (Wn1*Wn1 + Wn2*Wn2 + 4*xi1*xi2*Wn1*Wn2)*(Ix*Iz)
    LHS3 = (2*Wn1*Wn2*(xi1*Wn2 + xi2*Wn1))*(Ix*Iz)
    LHS4 = (Wn1*Wn1*Wn2*Wn2)*(Ix*Iz)
    diff1 = abs(RHS1 - LHS1)
    diff2 = abs(RHS2 - LHS2)
    diff3 = abs(RHS3 - LHS3)
    diff4 = abs(RHS4 - LHS4)

    print("RHS1 : ", RHS1, " LHS1 : ", LHS1, " Diff1 : ", diff1)
    print("RHS2 : ", RHS2, " LHS2 : ", LHS2, " Diff2 : ", diff2)
    print("RHS3 : ", RHS3, " LHS3 : ", LHS3, " Diff3 : ", diff3)
    print("RHS4 : ", RHS4, " LHS4 : ", LHS4, " Diff4 : ", diff4)


    if diff1 < e and diff2 < e and diff3 < e and diff4 < e :
        print("All equations satisfied")
    else :
        print("All equations are NOT satisfied")


def calculate_parameters(a):

    (Ix, Iy, Iz) = system_variables()

    wo = 2*pi/86400         # orbit frequency in rad/sec
    #wo = 7.236e-5
    xi1 = 0.05               # damping coefficient of closed loop nutation frequency poles
    xi2 = 0.7               # damping coefficient of closed loop orbit rate poles

    Tdx_max = 5e-6          # maximum magnitude of disturbance torque
    Tdz_max = 5e-6          # maximum magnitude of disturbance torque
    phi_ss = 0.05          # steady state error in roll in Deg
    psi_ss = 0.4           # steady state error in yaw in Deg

    phi_ss = radians(phi_ss)
    psi_ss = radians(psi_ss)

    # calculating h, kx without the approximation kx >> wo*h
    #Kx = (Tdx_max*(psi_ss/phi_ss) - Tdz_max)/(psi_ss - a*phi_ss)
    #h = (Tdx_max/phi_ss - Kx)/wo

    # calculating h, kx with approximation kx >> wo*h
    h = (Tdz_max + a*Tdx_max)/wo*psi_ss
    Kx = Tdx_max/phi_ss - wo*h



    print(h)
    print(Kx)


    P = wo*h*(Ix + Iz) + Iz*Kx + h*h
    Q = sqrt((wo*wo*h*h + wo*h*Kx)/(Ix*Iz))

    p = np.zeros(5)

    # defining p such that f(x) : p[0]*x**4 + p[1]*x**3 + p[2]*x**2 + p[3]*x + p[4].
    p[0] = 4*Ix*Ix*xi1*xi2*Q*Iz - 4*Ix*Ix*xi1*xi1*wo*h - Ix*Iz
    p[1] = 0
    p[2] = 4*Ix*Ix*xi1*xi1*Q*Q*Iz + 4*Ix*Ix*xi2*xi2*Q*Q*Iz - 8*Ix*Ix*xi1*wo*h*xi2*Q + P - 4*xi1*xi2*Ix*Iz*Q
    p[3] = 0
    p[4] = 4*Ix*Ix*xi1*xi2*Iz*Q*Q*Q - 4*Ix*Ix*xi2*xi2*Q*Q*wo*h - Ix*Iz*Q*Q

    quartic_roots = np.roots(p)

    print(quartic_roots)

    positive_roots = [x for x in quartic_roots if x > 0]

    if len(positive_roots) != 2:
        print("Error in closed loop frequency calculation")
        exit
    else : 
        Wn1 = max(positive_roots)
        Wn2 = min(positive_roots)
        Kxd = 2*Ix*(xi1*Wn1 + xi2*Wn2)
        a_psi = 2*Wn1*Wn2*(xi1*Wn2 + xi2*Wn1)*Ix*Iz - wo*h*Kxd


    calculated_parameters = (wo, xi1, xi2, h, Kx, Kxd, a_psi, Wn1, Wn2)

    return calculated_parameters


# set error tolerence in alpha
alpha_tolerance_deg = 1e-10            # tolerance in alpha in Deg
alpha_tolerance = radians(alpha_tolerance_deg)

alpha_guess_deg = 5
alpha_guess = radians(alpha_guess_deg)
a_guess = tan(alpha_guess)

simulation_parameters = calculate_parameters(a_guess)
a_calculated = simulation_parameters[6]
alpha_calculated = atan(a_calculated)

alpha_error = abs(alpha_calculated - alpha_guess)

while alpha_error > alpha_tolerance:

    alpha_guess = alpha_guess - (alpha_guess - alpha_calculated)/2
    a_guess = tan(alpha_guess)

    simulation_parameters = calculate_parameters(a_guess)
    print(len(simulation_parameters))
    a_calculated = simulation_parameters[6]
    alpha_calculated = atan(a_calculated)

    alpha_error = abs(alpha_calculated - alpha_guess)

display_parameters(simulation_parameters)
