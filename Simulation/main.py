from run_simulation import *
import math

# initializing simulation runtime parameters
start_time = 0                      # simulation start time in sec
end_time = 600                     # simulation end time in sec
time_step = 0.05                       # timestep Delta t
runtime_parameters = [start_time, end_time, time_step]

# initializing system variables
MOI_x = 1475.8425                   # moment of inertia about x-axis in kg.m^2
MOI_y = 886.2425                    # moment of inertia about x-axis in kg.m^2
MOI_z = 1532.96                     # moment of inertia about x-axis in kg.m^2

wheel_momentum = 20                 # angular momentum of momentum wheel in N.m.sec
orbit_rate = 2*math.pi/86164        # orbit rate in rad/sec

nutation_frequency = wheel_momentum/sqrt(MOI_x*MOI_z)
system_variables = [MOI_x, MOI_y, MOI_z, orbit_rate, nutation_frequency, wheel_momentum]

# initializing initial value of state variables
initial_roll = radians(0)     # initial roll angle in radians
initial_yaw = radians(0)        # initial yaw angle in radians
initial_roll_rate = 0           # initial roll rate in radians/sec
initial_yaw_rate = 0            # initial yaw rate in radians/sec
initial_conditions = [initial_roll, initial_roll_rate, initial_yaw, initial_yaw_rate]

# initializing System Control variables
roll_desired = radians(0.25)       # desired roll angle in radians (5 deg = 0.0872665)
offset_angle = 42             # offset angle of thruster in degree

# time between each measurement of earth sensor and control command
output_time = 0.5
controller_time = 0.5

sampling_parameters = [output_time, controller_time]

run_simulation(roll_desired, offset_angle, runtime_parameters, system_variables,sampling_parameters, initial_conditions)
