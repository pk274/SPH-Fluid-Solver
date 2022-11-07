// Paul Kull, 2022

#include <iostream>
#include <SFML/Graphics.hpp>

#include "./Simulation.h"

int main()
{
    // TestManager::test_kernel_integral();
    Simulation simulation = Simulation(SmallBox);
    simulation.run();
    return 0;
}