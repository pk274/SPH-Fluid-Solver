# Paul Kull, 2022

import matplotlib.pyplot as plt 
import numpy as np


timestep = [0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.0056, 0.006]
totaltime = [8.029, 6.725, 5.953, 5.648, 5.465, 5.340, 5.375, 5.627, 5.825, 5.943]
physicstime = [2.417, 2.357, 2.518, 2.750, 3.019, 3.204, 3.478, 3.891, 4.243, 4.478]
neighborhoodtime = [5.675, 4.368, 3.435, 2.898, 2.446, 2.136, 1.897, 1.736, 1.582, 1.465]



plt.plot(timestep, totaltime, label = 'Total time')
plt.plot(timestep, totaltime, 'bo')
plt.plot(timestep, physicstime, label = 'Pressure computation time')
plt.plot(timestep, physicstime, 'ro')
plt.plot(timestep, neighborhoodtime, label = 'Neighborhood search time')
plt.plot(timestep, neighborhoodtime, 'go')
plt.legend()
plt.title('Computation time for varying time step sizes')
plt.xlabel('time step size [s]')
plt.ylabel('computation time [s]')
plt.show()

