# Paul Kull, 2022

import matplotlib.pyplot as plt 
import matplotlib.ticker as plticker
import numpy as np


# Options:
plotAvgDensityAndIterations = 0
plotTimeStep = 0
plotDensityOnly = 1
plotDensityAndEstimatedDensity = 0
plotDensityAndAverageDensity = 0
plotCflNumber = 0

DensityThreshhold = 0.001
FrameIntervals = 1 / 25


# Plot AvgDensity and NumIterations
if (plotAvgDensityAndIterations):
    fig, (density, iters) = plt.subplots(2, sharex = True)
    avgDensityData = np.loadtxt("./avgDensityFile.dat", dtype=float)
    iterationsData = np.loadtxt("./iterationsFile.dat", dtype=float)

    iterationsSum = 0
    iterationsMax = 0
    for i in range(len(iterationsData)):
        iterationsSum += iterationsData[i][1]
        if iterationsData[i][1] > iterationsMax:
            iterationsMax = iterationsData[i][1]
    iterationsSum = iterationsSum / len(iterationsData)
    print('average number of iterations = ', iterationsSum)
    print('max number of iterations = ', iterationsMax)

    densitySum = 0
    for i in range(len(avgDensityData)):
        densitySum += avgDensityData[i][1]
    densitySum = densitySum / len(avgDensityData)
    print('mean average density = ', densitySum)

    iters.plot(*zip(*iterationsData), 'b-', label = '# iterations')
    #iters.plot([0, iterationsData[-1][0]], [iterationsSum, iterationsSum], 'm--', label = 'Average # iterations')

    density.plot(*zip(*avgDensityData))
    density.plot([0, avgDensityData[-1][0]], [1 + DensityThreshhold, 1 + DensityThreshhold], label = 'Density Threshhold')
    #fig.legend()
    iters.set_xlabel('simulated time (s)')
    iters.set_ylabel('# iterations')

    density.set_xlabel('simulated time (s)')
    density.set_ylabel('average density')
    plt.show()



if (plotTimeStep):
    timeStepData = np.loadtxt("./timeStepFile.dat", dtype=float)
    timeStepSum = 0
    for i in range(len(timeStepData)):
        timeStepSum += timeStepData[i][1]
    timeStepSum = timeStepSum / len(timeStepData)

    fig = plt.figure()
    timeStep = fig.add_subplot()

    #timeStep.plot(*zip(*timeStepData), 'g.', label = 'time step size')
    timeStep.plot(*zip(*timeStepData), 'g-', label = 'time step size')
    #timeStep.plot([0, timeStepData[-1][0]], [timeStepSum, timeStepSum], 'y--', label = 'Average time step size')
    #timeStep.set_ylim(0.003, 0.005)
    timeStep.set_title('Time step size')
    timeStep.set_xlabel('t')
    timeStep.set_ylabel('delta t')
    #Add Frame Grid
    loc = plticker.MultipleLocator(base=FrameIntervals)
    timeStep.yaxis.set_minor_locator(loc)
    timeStep.xaxis.set_minor_locator(loc)
    #timeStep.grid(which='minor', axis='both', linestyle='-', color = 'black', linewidth = '0.5')
    plt.show()



if (plotDensityOnly or plotDensityAndEstimatedDensity or plotDensityAndAverageDensity):
    avgDensityData = np.loadtxt("./avgDensityFile.dat", dtype=float)
    densitySum = 0
    for i in range(len(avgDensityData)):
        densitySum += avgDensityData[i][1]
    densitySum = densitySum / len(avgDensityData)
    plt.plot(*zip(*avgDensityData))
    if (plotDensityAndAverageDensity):
        plt.plot([0, avgDensityData[-1][0]], [densitySum, densitySum], 'g--', label = 'mean average density')
    #else:
    #    plt.plot([0, avgDensityData[-1][0]], [1 + DensityThreshhold, 1 + DensityThreshhold], label = 'density threshhold')

    if (plotDensityAndEstimatedDensity):
        estimatedDensityData = np.loadtxt("./estimatedDensityFile.dat", dtype=float)
        plt.plot(*zip(*estimatedDensityData), 'g-', label = 'estimated average density')
    #plt.legend(["Average density", "Density threshhold", "Estimated average density"])
    plt.xlabel('simulated time (s)')
    plt.ylabel('average density')
    plt.show()


if (plotCflNumber):
    cflNumberData = np.loadtxt("./cflNumberFile.dat", dtype=float)
    cflNumberSum = 0
    for i in range(len(cflNumberData)):
        cflNumberSum += cflNumberData[i][1]
    cflNumberSum = cflNumberSum / len(cflNumberData)
    plt.plot(*zip(*cflNumberData), 'm-')
    plt.plot([0, cflNumberData[-1][0]], [cflNumberSum, cflNumberSum], 'y--', label = 'average CFL number')
    plt.legend(["CFL number", "Average CFL number"])
    plt.xlabel('simulated time (s)')
    plt.ylabel('CFL number')
    plt.show()
