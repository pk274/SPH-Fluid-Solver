# Paul Kull, 2022

import matplotlib.pyplot as plt 
import matplotlib.ticker as plticker
import numpy as np

timeStepData = np.loadtxt("./timeStepFile.dat", dtype=float)
iterationsData = np.loadtxt("./iterationsFile.dat", dtype=float)

timeStepSum = 0
for i in range(len(timeStepData)):
    timeStepSum += timeStepData[i][1]
timeStepSum = timeStepSum / len(timeStepData)

iterationsSum = 0
for i in range(len(iterationsData)):
    iterationsSum += iterationsData[i][1]
iterationsSum = iterationsSum / len(iterationsData)



fig, (timeStep, iters) = plt.subplots(2)
timeStep.plot(*zip(*timeStepData), 'g.', label = 'time step size')
timeStep.plot(*zip(*timeStepData), 'g-', label = 'time step size')
timeStep.plot([0, timeStepData[-1][0]], [timeStepSum, timeStepSum], 'y--', label = 'Average time step size')
#timeStep.set_ylim(0.003, 0.005)
iters.plot(*zip(*iterationsData), 'b.', label = '# iterations')
iters.plot(*zip(*iterationsData), 'b-', label = '# iterations')
iters.plot([0, iterationsData[-1][0]], [iterationsSum, iterationsSum], 'm--', label = 'Average # iterations')
#fig.legend()
timeStep.set_title('Time step size')
iters.set_title('# iterations')
timeStep.set_xlabel('t')
timeStep.set_ylabel('delta t')
#Add Frame Grid
intervals = 1 / 30
loc = plticker.MultipleLocator(base=intervals)
timeStep.yaxis.set_minor_locator(loc)
timeStep.xaxis.set_minor_locator(loc)
timeStep.grid(which='minor', axis='both', linestyle='-', color = 'black', linewidth = '0.5')

iters.set_xlabel('t')
iters.set_ylabel('# iterations')
plt.show()


