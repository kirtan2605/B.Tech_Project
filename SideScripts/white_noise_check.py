import math
import numpy as np
import matplotlib.pyplot as plt
from pandas import Series
from pandas.plotting import autocorrelation_plot

T_sam = 0.1
run_time = 100
num_of_samples = int(run_time / T_sam) + 1
sampling_freq = 1 / T_sam
error_std_dev_deg = 0.012
error_std_dev_rad = math.radians(error_std_dev_deg)
white_noise = [np.random.normal(0, error_std_dev_rad, 1) for i in range(num_of_samples)]
White_Noise = Series(white_noise)
x_axis = range(num_of_samples)

plt.plot(x_axis, white_noise)
plt.title("White Noise (\u03C3 = 0.012 Deg = 0.00021 rad)", fontsize=12)
plt.xlabel("Sample", fontsize=12)
plt.ylabel("Noise (rad)", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.show()

# check with Autocorrelation Plots
autocorrelation_plot(White_Noise)
plt.show()

## check with Power Spectral Density
#plt.psd(White_Noise, Fs=sampling_freq)
#plt.xlim([0, 20])
#plt.show()
