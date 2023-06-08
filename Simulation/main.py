from run_simulation import *
import math

# initializing simulation runtime parameters
start_time = 0                          # simulation start time in sec
end_time = 600                         # simulation end time in sec
time_step = 0.1                         # timestep Delta t
runtime_parameters = [start_time, end_time, time_step]

# initializing initial value of state variables
initial_roll = radians(0)               # initial roll angle in radians
initial_yaw = radians(0)                # initial yaw angle in radians
initial_roll_rate = radians(0)          # initial roll rate in radians/sec
initial_yaw_rate = radians(0)           # initial yaw rate in radians/sec
initial_conditions = [initial_roll, initial_roll_rate, initial_yaw, initial_yaw_rate]

# initializing System Control variables
roll_desired = radians(5)       # desired roll angle in radians (5 deg = 0.0872665)
offset_angle = 41.8                # offset angle of thruster in degree

# time between each measurement of earth sensor and control command
output_time = 0.5
controller_time = 0.5

sampling_parameters = [output_time, controller_time]

run_simulation(roll_desired, offset_angle, runtime_parameters, sampling_parameters, initial_conditions)
