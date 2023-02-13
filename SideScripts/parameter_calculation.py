from math import *


def calculated_a(a_assumed):
    ## spacecraft moments of inertia
    # calculating (INSAT-3DR : [2211, 2.4, 1.6, 1.5])
    mass = 2211             # mass in kg
    height = 2.4            # height in m
    width = 1.6             # width in m
    depth = 1.5             # depth in m
    # by convention height > width > depth
    Iz = mass*(height**2 + width**2)/12
    Ix = mass*(height**2 + depth**2)/12
    Iy = mass*(width**2 + depth**2)/12

    ## or define the moments of inertia
    # Ix = 800
    # Iy = 680
    # Iz = 1000

    print()
    print("Ix : ", Ix, " kg-m^2")
    print("Iy : ", Iy, " kg-m^2")
    print("Iz : ", Iz, " kg-m^2")
    print()

    wo = 2*pi/86400     # orbit frequency in rad/sec
    ze = 0.9           # damping coefficient of closed loop poles

    print("wo : ", wo)

    Tdx_max = 3e-6           # maximum magnitude of disturbance torque
    Tdz_max = 3e-6           # maximum magnitude of disturbance torque
    phi_ss = 0.045          # steady state error in roll in Deg
    phi_ss = radians(phi_ss)
    psi_ss = 0.25           # steady state error in yaw in Deg
    psi_ss = radians(psi_ss)


    # calculating h, kx without the approximation kx>>wo*h
    kx = (Tdx_max*(psi_ss/phi_ss) - Tdz_max)/(psi_ss - a_assumed*phi_ss)
    h = (Tdx_max/phi_ss - kx)/wo


    # calculating kxd,a_calculated, wn1, wn2 (THIS PART IS PERFECT)
    A = sqrt(((wo * wo * h * h) + (wo * h * kx)) / (Ix * Iz))
    B = 1 / (2 * Ix * ze)
    C = (wo * h * (Ix + Iz) + Iz * kx + h * h) / (Ix * Iz)

    Kxd = sqrt(((C + A * (2 - 4 * ze * ze)) * kx) / ((kx * B * B) - (2 * ze * A * B) + (h * wo) / (Ix * Iz)))
    a_psi = (2 * ze * A * B * Kxd * Ix * Iz - wo * Kxd * h) / (h * kx)
    Wn1 = ((Kxd * B) + sqrt(Kxd * Kxd * B * B - 4 * A)) / 2
    Wn2 = A / Wn1

    print("h :", h)
    print("Kx: " + str(kx))
    print("Kxd: " + str(Kxd))
    print("a_psi: " + str(a_psi))
    print("Wn1: " + str(Wn1))
    print("Wn2: " + str(Wn2))
    print()



    return a_psi

# the initial assumption of the value of a starts with 1
# corrections are made in the value of a till a desired accuracy is achieved


# set error tolerence in alpha
alpha_tolerance_deg = 0.0005          # tolerance in alpha in Deg
alpha_tolerance = alpha_tolerance_deg*pi/180

alpha_guess_deg = 5
alpha_guess = radians(alpha_guess_deg)
a_guess = tan(alpha_guess)


a_return = calculated_a(a_guess)
alpha_calculated = atan(a_return)

alpha_error = abs(alpha_calculated - alpha_guess)

while alpha_error > alpha_tolerance:

    alpha_guess = alpha_guess - (alpha_guess - alpha_calculated)/10

    a_guess = tan(alpha_guess)

    a_return = calculated_a(a_guess)
    alpha_calculated = atan(a_return)

    alpha_error = abs(alpha_calculated - alpha_guess)

print("a guess : ", a_guess)
print("a returned : ", a_return)



