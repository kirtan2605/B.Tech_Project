import numpy as np
import math

"""
    Adapted from : 
    https://medium.com/geekculture/runge-kutta-numerical-integration-of-ordinary-differential-equations-in-python-9c8ab7fb279c
"""

def ode_system(_t, _x):
    """
    system of first order differential equations
    _t: discrete time step value
    _x: state vector [x1, x2]
    """
    # initializing system variables
    Ix = 1475.8425                   # moment of inertia about x-axis in kg.m^2
    Iy = 886.2425                    # moment of inertia about x-axis in kg.m^2
    Iz = 1532.96                     # moment of inertia about x-axis in kg.m^2

    h = 20                 # angular momentum of momentum wheel in N.m.sec
    omega0 = 2*math.pi/86164        # orbit rate in rad/sec

    a = 4 * (omega0 ** 2) * (Iy - Iz)
    b = -1 * omega0 * (Ix - Iy + Iz)
    c = (omega0 ** 2) * (Iy - Ix)

    Tx = 0
    Tz = 0

    return np.array([_x[1], (Tx - (a + omega0*h)*_x[0] - (b + h)*_x[3])/Ix, _x[3], (Tz - (c + omega0*h)*_x[2] + (b + h)*_x[1])/Ix ])


# runge-kutta fourth-order numerical integration
def rk4(func, tk, _yk, _dt=0.01, **kwargs):
    """
    single-step fourth-order numerical integration (RK4) method
    func: system of first order ODEs
    tk: current time step
    _yk: current state vector [y1, y2, y3, ...]
    _dt: discrete time step size
    **kwargs: additional parameters for ODE system
    returns: y evaluated at time k+1
    """

    # evaluate derivative at several stages within time interval
    f1 = func(tk, _yk, **kwargs)
    f2 = func(tk + _dt / 2, _yk + (f1 * (_dt / 2)), **kwargs)
    f3 = func(tk + _dt / 2, _yk + (f2 * (_dt / 2)), **kwargs)
    f4 = func(tk + _dt, _yk + (f3 * _dt), **kwargs)

    # return an average of the derivative over tk, tk + dt
    return _yk + (_dt / 6) * (f1 + (2 * f2) + (2 * f3) + f4)