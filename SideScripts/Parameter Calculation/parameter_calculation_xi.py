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
    xi = 0.7                # damping coefficient of closed loop poles

    Tdx_max = 5e-6          # maximum magnitude of disturbance torque
    Tdz_max = 5e-6          # maximum magnitude of disturbance torque
    phi_ss = 0.05          # steady state error in roll in Deg
    psi_ss = 0.4           # steady state error in yaw in Deg

    phi_ss = radians(phi_ss)
    psi_ss = radians(psi_ss)


    # calculating h, kx without the approximation kx >> wo*h !!
    #Kx = (Tdx_max*(psi_ss/phi_ss) - Tdz_max)/(psi_ss - a_assumed*phi_ss)
    #h = (Tdx_max/phi_ss - Kx)/wo

    # calculating h, kx with the approximation kx >> wo*h
    h = (Tdz_max - a_assumed*Tdx_max)/(wo*psi_ss)
    Kx = Tdx_max/phi_ss - wo*h

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
    Ix,Iy,Iz = system_variables()
    print("\nSystem Parameters")
    print("Ix :", Ix)
    print("Iy :", Iy)
    print("Iz :", Iz , "\n")

    wo, xi, h, Kx, Kxd, a_psi, Wn1, Wn2 = parameters
    print("\nSimulation Parameters")
    print("w0 :", wo)
    print("xi :", xi)
    print("h :", h)
    print("Kx: ", Kx)
    print("Kxd: ", Kxd)
    print("a_psi: ", a_psi)
    print("alpha: ", degrees(atan(a_psi)))
    print("Wn1: ", Wn1)
    print("Wn2: ", Wn2 , "\n")      

    # verify if the equations hold for epsilon e
    e = 1
    
    RHS1 = (Iz*Kxd)
    RHS2 = (wo*h*(Ix+Iz) + Iz*Kx + h*h + a_psi*h*Kxd)
    RHS3 = (a_psi*h*Kx + wo*h*Kxd)
    RHS4 = (wo*wo*h*h + wo*h*Kx)
    LHS1 = (2*(xi*Wn1 + xi*Wn2))*(Ix*Iz)
    LHS2 = (Wn1*Wn1 + Wn2*Wn2 + 4*xi*xi*Wn1*Wn2)*((Ix*Iz))
    LHS3 = (2*Wn1*Wn2*xi*(Wn1 + Wn2) )*(Ix*Iz)
    LHS4 = (Wn1*Wn1*Wn2*Wn2)*(Ix*Iz)
    diff1 = abs(RHS1 - LHS1)
    diff2 = abs(RHS2 - LHS2)
    diff3 = abs(RHS3 - LHS3)
    diff4 = abs(RHS4 - LHS4)

    diff1RHSp = (diff1/RHS1)*100
    diff1LHSp = (diff1/LHS1)*100
    diff2RHSp = (diff2/RHS2)*100
    diff2LHSp = (diff2/LHS2)*100
    diff3RHSp = (diff3/RHS3)*100
    diff3LHSp = (diff3/LHS3)*100
    diff4RHSp = (diff4/RHS4)*100
    diff4LHSp = (diff4/LHS4)*100

    print("RHS1 percentage : ", diff1RHSp,"%")
    print("LHS1 percentage : ", diff1LHSp,"%")
    if diff1RHSp < e and diff1LHSp < e :
        print("Equation 1 satisfied")
    else :
        print("Equation 1 NOTe satisfied")
    print()
    
    print("RHS2 percentage : ", diff2RHSp,"%")
    print("LHS2 percentage : ", diff2LHSp,"%")
    if diff2RHSp < e and diff2LHSp < e :
        print("Equation 2 satisfied")
    else :
        print("Equation 2 NOT satisfied")
    print()

    print("RHS3 percentage : ", diff3RHSp,"%")
    print("LHS3 percentage : ", diff3LHSp,"%")
    if diff3RHSp < e and diff3LHSp < e :
        print("Equation 3 satisfied")
    else :
        print("Equation 3 NOT satisfied")
    print()

    print("RHS4 percentage : ", diff4RHSp,"%")
    print("LHS4 percentage : ", diff4LHSp,"%")
    if diff4RHSp < e and diff4LHSp < e :
        print("Equation 4 satisfied")
    else :
        print("Equation 4 NOT satisfied")
    print()

# set error tolerence in alpha
alpha_tolerance_deg = 1e-10            # tolerance in alpha in Deg
alpha_tolerance = radians(alpha_tolerance_deg)

alpha_guess_deg = 5
alpha_guess = radians(alpha_guess_deg)
a_guess = tan(alpha_guess)

simulation_parameters = calculate_parameters(a_guess)
a_calculated = simulation_parameters[5]
alpha_calculated = atan(a_calculated)

alpha_error = abs(alpha_calculated - alpha_guess)

while alpha_error > alpha_tolerance:

    alpha_guess = alpha_guess - (alpha_guess - alpha_calculated)/2
    a_guess = tan(alpha_guess)

    simulation_parameters = calculate_parameters(a_guess)
    a_calculated = simulation_parameters[5]
    alpha_calculated = atan(a_calculated)

    alpha_error = abs(alpha_calculated - alpha_guess)

display_parameters(simulation_parameters)