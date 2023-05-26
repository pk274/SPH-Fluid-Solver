# Paul Kull, 2022

import matplotlib.pyplot as plt 
import numpy as np

#close particles
#particles = [512, 732, 952, 1174]
#uniformGrid = [10.8595, 19.9138, 27.4366, 35.6074732]
#improvedGrid = [10.617, 18.951, 28.5852, 36.171]
#nSquare = [9.933, 22.011, 31.7204, 57.768]

#bd small
#particles = [412, 712, 1041, 1212]
#uniformGrid = [7.046, 19.093, 30.697, 38.142]
#improvedGrid = [6.312, 17.971, 29.1528, 36.191]
#nSquare = [5.753, 20.095, 37.3574, 60.471]

#bd big
particles = [732, 1032, 1532, 2232, 3132]
particlesNSquare = [732, 1032, 1532, 2232]
uniformGrid = [14.693, 27.905, 48.593, 78.599, 118.5342]
improvedGrid = [10.098, 22.679, 43.675, 71.765, 108.999]
nSquare = [10.443, 28.623, 75.543, 174.819]



plt.plot(particlesNSquare, nSquare, 'ro', label = 'naive search')
plt.plot(particles, uniformGrid, 'bo', label = 'uniform grid')
plt.plot(particles, improvedGrid, 'go', label = 'improved grid')
plt.legend()
plt.title('Time spent with neighborhood search')
plt.xlabel('# of particles')
plt.ylabel('time in ms')
plt.show()

