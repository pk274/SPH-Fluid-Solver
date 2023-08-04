# Paul Kull, 2022

import matplotlib.pyplot as plt 
import numpy as np


timeStep = [0.0045, 0.004, 0.003, 0.0025, 0.002, 0.0015, 0.001, 0.0008, 0.0006, 0.0004]
maxAvgError = [0.132, 0.099, 0.059, 0.0425, 0.029, 0.0175, 0.00825, 0.00565, 0.00305, 0.00149]
meanAvgError = [0.051, 0.04, 0.023, 0.0158, 0.01, 0.0055, 0.0025, 0.00165, 0.00095, 0.0004]
timeStep.reverse()
maxAvgError.reverse()
meanAvgError.reverse()


plt.plot(timeStep, maxAvgError, 'r-', label = 'maxAvgError')
plt.plot(timeStep, meanAvgError, 'c-', label = 'meanAvgError')
plt.legend(['max avg. density error', 'mean avg. density error'])
plt.xlabel('time step size')
plt.ylabel('density error')
plt.show()

