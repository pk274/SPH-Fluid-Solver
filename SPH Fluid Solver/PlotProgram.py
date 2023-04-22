# Paul Kull, 2022

import matplotlib.pyplot as plt 
import numpy as np

avgDensityData = np.loadtxt("./avgDensityFile.dat", dtype=float)
sum = 0
for i in range(len(avgDensityData)):
    sum += avgDensityData[i][1]
sum = sum / len(avgDensityData)
plt.plot(*zip(*avgDensityData))
#plt.plot([0, avgDensityData[-1][0]], [sum, sum], label = 'Mean average density')
plt.plot([0, avgDensityData[-1][0]], [1.01, 1.01], label = 'Density Threshhold')#
plt.legend()
plt.title('Average density over time')
plt.xlabel('t')
plt.ylabel('Average Density')
plt.show()

