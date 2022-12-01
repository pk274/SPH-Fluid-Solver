# Paul Kull, 2022

import matplotlib.pyplot as plt 
import numpy as np

avgDensityData = np.loadtxt("./avgDensityFile", dtype=float)
plt.plot(*zip(*avgDensityData))
plt.title('Average density over time')
plt.xlabel('t')
plt.ylabel('p')
plt.legend(['Average density'])
plt.show()
