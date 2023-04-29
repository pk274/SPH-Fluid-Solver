# Paul Kull, 2022

import matplotlib.pyplot as plt 
import matplotlib.ticker as plticker
import numpy as np


# Options:
plotAvgDensityAndIterations = True
plotTimeStep = True
plotDensityOnly = False

DensityThreshhold = 0.01
FrameIntervals = 1 / 25
# Plot AvgDensity and NumIterations


if (plotAvgDensityAndIterations):
    fig, (density, iters) = plt.subplots(2, sharex = True)
    avgDensityData = np.loadtxt("./avgDensityFile.dat", dtype=float)
    iterationsData = np.loadtxt("./iterationsFile.dat", dtype=float)

    iterationsSum = 0
    for i in range(len(iterationsData)):
        iterationsSum += iterationsData[i][1]
    iterationsSum = iterationsSum / len(iterationsData)

    densitySum = 0
    for i in range(len(avgDensityData)):
        densitySum += avgDensityData[i][1]
    densitySum = densitySum / len(avgDensityData)

    iters.plot(*zip(*iterationsData), 'b-', label = '# iterations')
    iters.plot([0, iterationsData[-1][0]], [iterationsSum, iterationsSum], 'm--', label = 'Average # iterations')

    density.plot(*zip(*avgDensityData))
    density.plot([0, avgDensityData[-1][0]], [1 + DensityThreshhold, 1 + DensityThreshhold], label = 'Density Threshhold')
    #fig.legend()
    iters.set_title('# iterations')
    iters.set_xlabel('t')
    iters.set_ylabel('# iterations')

    density.set_title('Average density over time')
    density.set_xlabel('t')
    density.set_ylabel('Average Density')
    plt.show()



if (plotTimeStep):
    timeStepData = np.loadtxt("./timeStepFile.dat", dtype=float)
    timeStepSum = 0
    for i in range(len(timeStepData)):
        timeStepSum += timeStepData[i][1]
    timeStepSum = timeStepSum / len(timeStepData)

    fig = plt.figure()
    timeStep = fig.add_subplot()

    timeStep.plot(*zip(*timeStepData), 'g.', label = 'time step size')
    timeStep.plot(*zip(*timeStepData), 'g-', label = 'time step size')
    timeStep.plot([0, timeStepData[-1][0]], [timeStepSum, timeStepSum], 'y--', label = 'Average time step size')
    #timeStep.set_ylim(0.003, 0.005)
    timeStep.set_title('Time step size')
    timeStep.set_xlabel('t')
    timeStep.set_ylabel('delta t')
    #Add Frame Grid
    loc = plticker.MultipleLocator(base=FrameIntervals)
    timeStep.yaxis.set_minor_locator(loc)
    timeStep.xaxis.set_minor_locator(loc)
    timeStep.grid(which='minor', axis='both', linestyle='-', color = 'black', linewidth = '0.5')
    plt.show()



if (plotDensityOnly):
    avgDensityData = np.loadtxt("./avgDensityFile.dat", dtype=float)
    densitySum = 0
    for i in range(len(avgDensityData)):
        densitySum += avgDensityData[i][1]
    densitySum = densitySum / len(avgDensityData)
    plt.plot(*zip(*avgDensityData))
    #plt.plot([0, avgDensityData[-1][0]], [sum, sum], label = 'Mean average density')
    plt.plot([0, avgDensityData[-1][0]], [DensityThreshhold, DensityThreshhold], label = 'Density Threshhold')
    plt.legend()
    plt.title('Average density over time')
    plt.xlabel('t')
    plt.ylabel('Average Density')
    plt.show()
