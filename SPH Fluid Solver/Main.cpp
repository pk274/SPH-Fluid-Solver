// Paul Kull, 2022

#include <iostream>
#include <SFML/Graphics.hpp>

#include "./Simulation.h"


int main()
{
    Simulation simulation = Simulation(RotatedBox);
    simulation.run();
    return 0;
}
