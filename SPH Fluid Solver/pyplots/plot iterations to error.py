# Paul Kull, 2022

import matplotlib.pyplot as plt 
import numpy as np


iterations = [2, 3, 4, 5, 7, 8, 12, 16, 24, 32]
maxAvgError = [0.019, 0.013, 0.0105, 0.0085, 0.0064, 0.0055, 0.0036, 0.0026, 0.00178, 0.00129]
meanAvgError = [0.0065, 0.0042, 0.0033, 0.0027, 0.00188, 0.0016, 0.00105, 0.0006, 0.00052, 0.00039]




plt.plot(iterations, maxAvgError, 'r-', label = 'maxAvgError')
plt.plot(iterations, meanAvgError, 'c-', label = 'meanAvgError')
plt.legend(['max avg. density error', 'mean avg. density error'])
plt.xlabel('# of iterations')
plt.ylabel('density error')
plt.show()

