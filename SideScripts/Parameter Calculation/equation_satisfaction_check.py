from math import *

def check_equation_satisfaction(a_psi, h, Kx, Kxd, Wn1, Wn2, xi1, xi2):
    Ix = 800
    Iz = 1000

    wo = 7.236e-5

    Tdx = 5e-6
    Tdz = 5e-6

    phi_ss = radians(0.05)
    psi_ss = radians(0.4)

    # verify if the equations hold for epsilon e
    e = 5
    
    RHS1 = phi_ss
    LHS1 = Tdx/(wo*h + Kx)

    RHS2 = psi_ss
    # No approximation LHS
    LHS2 = (Tdz/(wo*h)) - (a_psi*Tdx*Kx)/(wo*h*(Kx + wo*h))
    # Approximation LHS
    #LHS2 = ((Tdz - a_psi*Tdx)/(wo*h))

    RHS3 = (Iz*Kxd)
    LHS3 = (2*(xi1*Wn1 + xi2*Wn2))*(Ix*Iz)

    RHS4 = (wo*h*(Ix+Iz) + Iz*Kx + h*h + a_psi*h*Kxd)
    LHS4 = (Wn1*Wn1 + Wn2*Wn2 + 4*xi1*xi2*Wn1*Wn2)*(Ix*Iz)

    RHS5 = (a_psi*h*Kx + wo*h*Kxd)
    LHS5 = (2*Wn1*Wn2*(xi1*Wn2 + xi2*Wn1))*(Ix*Iz)

    RHS6 = (wo*wo*h*h + wo*h*Kx)
    LHS6 = (Wn1*Wn1*Wn2*Wn2)*(Ix*Iz)

    diff1 = abs(RHS1 - LHS1)
    diff2 = abs(RHS2 - LHS2)
    diff3 = abs(RHS3 - LHS3)
    diff4 = abs(RHS4 - LHS4)
    diff3 = abs(RHS5 - LHS5)
    diff4 = abs(RHS6 - LHS6)

    diff1RHSp = (diff1/RHS1)*100
    diff1LHSp = (diff1/LHS1)*100
    diff2RHSp = (diff2/RHS2)*100
    diff2LHSp = (diff2/LHS2)*100
    diff3RHSp = (diff3/RHS3)*100
    diff3LHSp = (diff3/LHS3)*100
    diff4RHSp = (diff4/RHS4)*100
    diff4LHSp = (diff4/LHS4)*100
    diff5RHSp = (diff3/RHS5)*100
    diff5LHSp = (diff3/LHS5)*100
    diff6RHSp = (diff4/RHS6)*100
    diff6LHSp = (diff4/LHS6)*100

    print("RHS1 percentage : ", diff1RHSp,"%")
    print("LHS1 percentage : ", diff1LHSp,"%")

    print("RHS2 percentage : ", diff2RHSp,"%")
    print("LHS2 percentage : ", diff2LHSp,"%")

    print("RHS3 percentage : ", diff3RHSp,"%")
    print("LHS3 percentage : ", diff3LHSp,"%")

    print("RHS4 percentage : ", diff4RHSp,"%")
    print("LHS4 percentage : ", diff4LHSp,"%")

    print("RHS5 percentage : ", diff5RHSp,"%")
    print("LHS5 percentage : ", diff5LHSp,"%")

    print("RHS6 percentage : ", diff6RHSp,"%")
    print("LHS6 percentage : ", diff6LHSp,"%")

    print()

    if diff1RHSp < e and diff1LHSp < e :
        print("Equation 1 satisfied")
    else :
        print("Equation 1 NOT satisfied")
    print()

    if diff2RHSp < e and diff2LHSp < e :
        print("Equation 2 satisfied")
    else :
        print("Equation 2 NOT satisfied")
    print()

    if diff3RHSp < e and diff3LHSp < e :
        print("Equation 3 satisfied")
    else :
        print("Equation 3 NOT satisfied")
    print()

    if diff4RHSp < e and diff4LHSp < e :
        print("Equation 4 satisfied")
    else :
        print("Equation 4 NOT satisfied")
    print()

    if diff5RHSp < e and diff5LHSp < e :
        print("Equation 5 satisfied")
    else :
        print("Equation 5 NOT satisfied")
    print()

    if diff6RHSp < e and diff6LHSp < e :
        print("Equation 6 satisfied")
    else :
        print("Equation 6 NOT satisfied")
    print()


# Case1 Check
check_equation_satisfaction(0.888, 20, 4.27e-3, 42.8, 8.48e-5, 0.0381, 0.7, 0.7)