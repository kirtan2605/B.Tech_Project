import math

# initializing system variables
Ix = 1475.8425                   # moment of inertia about x-axis in kg.m^2
Iy = 886.2425                    # moment of inertia about x-axis in kg.m^2
Iz = 1532.96                     # moment of inertia about x-axis in kg.m^2

h = 20                 # angular momentum of momentum wheel in N.m.sec
omega0 = 2*math.pi/86164        # orbit rate in rad/sec

a = 4 * (omega0 ** 2) * (Iy - Iz)
b = -1 * omega0 * (Ix - Iy + Iz)
c = (omega0 ** 2) * (Iy - Ix)

print('h : ', h)
print('a/wo : ', math.fabs(a/omega0))
print('b : ', math.fabs(b))
print('c/wo : ', math.fabs(c/omega0))
print('wo*(Ix + Iz) : ', math.fabs(omega0*(Ix+Iz)))