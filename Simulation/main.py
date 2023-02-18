from run_simulation import *
import math

# initializing simulation runtime parameters
start_time = 0              # simulation start time in sec
end_time = 5000             # simulation end time in sec
time_step = 0.5              # timestep Delta t
# time_step is the sampling time for our simulation
runtime_parameters = [start_time, end_time, time_step]

# initializing system variables
MOI_x = 1475.8425                   # moment of inertia about x-axis in kg.m^2
MOI_y = 886.2425                    # moment of inertia about x-axis in kg.m^2
MOI_z = 1532.96                     # moment of inertia about x-axis in kg.m^2
wheel_momentum = 11                  # angular momentum of momentum wheel in N.m.sec
orbit_rate = 2*math.pi/86400        # orbit rate in rad/sec

nutation_frequency = wheel_momentum/((MOI_x*MOI_z)**0.5)
system_variables = [MOI_x, MOI_y, MOI_z, orbit_rate, nutation_frequency, wheel_momentum]

# considering the system to be under

# initializing initial value of state variables
initial_roll = 0.0174533        # initial roll angle in radians
initial_yaw = 0                 # initial yaw angle in radians
initial_roll_rate = 0           # initial roll rate in radians/sec
initial_yaw_rate = 0            # initial yaw rate in radians/sec
initial_conditions = [initial_roll, initial_roll_rate, initial_yaw, initial_yaw_rate]

# initial condition on torques are specified in run_simulation.py

# initializing System Control variables
roll_desired = 0          # desired roll angle in radians (5 deg = 0.0872665)
offset_angle = 35.25          # offset angle of thruster in degree

run_simulation(roll_desired, offset_angle, runtime_parameters, system_variables, initial_conditions)
