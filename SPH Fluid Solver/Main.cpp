// Paul Kull, 2022

#include <iostream>
#include <SFML/Graphics.hpp>

#include "./Simulation.h"


int main()
{
    Simulation simulation = Simulation(ManyLayers);
    simulation.run();
    // simulation.render_from_file("testrender");
    return 0;
}
