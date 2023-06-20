# Paul Kull, 2022

import matplotlib.pyplot as plt 
import numpy as np

avgDensityDataIISPH = np.loadtxt("./avgDensityFileIISPH.txt", dtype=float)
avgDensityDataSESPH = np.loadtxt("./avgDensityFileSESPH.txt", dtype=float)

DensityThreshhold = 0.01


densitySumIISPH = 0
for i in range(len(avgDensityDataIISPH)):
    densitySumIISPH += avgDensityDataIISPH[i][1]
densitySumIISPH = densitySumIISPH / len(avgDensityDataIISPH)
plt.plot(*zip(*avgDensityDataIISPH))
densitySumSESPH = 0
for i in range(len(avgDensityDataSESPH)):
    densitySumSESPH += avgDensityDataSESPH[i][1]
densitySumSESPH = densitySumSESPH / len(avgDensityDataSESPH)
plt.plot(*zip(*avgDensityDataSESPH))
#plt.plot([0, avgDensityDataIISPH[-1][0]], [1 + DensityThreshhold, 1 + DensityThreshhold], label = 'density threshhold')

plt.legend(["Avg. density ISPH", "Avg. density SESPH", "Density threshold"])
plt.xlabel('simulated time (s)')
plt.ylabel('average density')
plt.show()
