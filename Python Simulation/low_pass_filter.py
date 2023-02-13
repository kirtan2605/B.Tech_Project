import numpy as np


def low_pass_filter(y, x, b, a, N):

    # adding zero pad when past data is insufficient
    if len(y) < (N + 1):
        y = np.pad(y, (N+1-len(y), 0), 'constant', constant_values=(0, 0))
    if len(x) < (N + 1):
        x = np.pad(x, (N+1-len(x), 0), 'constant', constant_values=(0, 0))

    low_pass_return_value = np.matmul(b, np.flip(x)) - np.matmul(a, np.flip(y))

    return low_pass_return_value
