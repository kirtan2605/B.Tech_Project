from run_simulation import *


# initializing simulation runtime parameters
start_time = 0              # simulation start time in sec
end_time = 500             # simulation end time in sec
time_step = 0.01              # timestep Delta t
# time_step is the sampling time for our simulation
# sampling time can be optimally set using Shannon Sampling Theorem
# to apply the sampling theorem, we need information about the bandwidth of the time response
# the time response can be analyzed from the root locus ??
runtime_parameters = [start_time, end_time, time_step]

# initializing system variables
MOI_x = 1355.818                    # moment of inertia about x-axis in kg.m^2
MOI_y = 1355.818                    # moment of inertia about x-axis in kg.m^2
MOI_z = 1355.818                    # moment of inertia about x-axis in kg.m^2
wheel_momentum = 203.2737           # angular momentum of momentum wheel in N.m.sec
orbit_rate = 0.000073               # orbit rate in rad/sec
nutation_frequency = wheel_momentum/((MOI_x*MOI_z)**0.5)
system_variables = [MOI_x, MOI_y, MOI_z, orbit_rate, wheel_momentum]

print("\nOrbit Rate : ", orbit_rate, "Hz")
print("Nutation Frequency : ", nutation_frequency, "Hz")
print("Nutation Frequency / Orbit Rate : ", nutation_frequency/orbit_rate)


print("\nSampling Frequency : ", 1/time_step, "Hz")
print("No information is lost due to sampling")
print("since f_max i.e nutation_rate < Sampling Frequency")


# considering the system to be under

# initializing initial value of state variables
initial_roll = 0.0174533        # initial roll angle in radians
initial_yaw = 0                 # initial yaw angle in radians
initial_roll_rate = 0           # initial roll rate in radians/sec
initial_yaw_rate = 0          # initial yaw rate in radians/sec
initial_conditions = [initial_roll, initial_roll_rate, initial_yaw, initial_yaw_rate]

# initial condition on torques are specified in run_simulation.py

# initializing System Control variables
roll_desired = 0          # desired roll angle in radians
offset_angle = 5          # offset angle of thruster in degree

run_simulation(roll_desired, offset_angle, runtime_parameters, system_variables, initial_conditions)
