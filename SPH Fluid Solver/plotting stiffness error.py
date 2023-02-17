# Paul Kull, 2022

import matplotlib.pyplot as plt 
import numpy as np

meanDensityError = [9.98, 4.1, 3, 2.21, 1.5, 1.2, 1, 0.685, 0.51, 0.345]
k = [100000, 250000, 325000, 500000, 650000, 850000, 1000000, 1500000, 2000000, 3000000]
plt.plot(k, meanDensityError)
plt.plot(k, meanDensityError, 'bo')
plt.title('Relation between k and the density error')
plt.xlabel('Stiffness k')
plt.ylabel('Mean average density error in percent')
plt.legend(['Density error in %'])
plt.show()
