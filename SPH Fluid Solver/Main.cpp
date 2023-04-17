// Paul Kull, 2022

#include <iostream>
#include <SFML/Graphics.hpp>

#include "./Simulation.h"
#include "./InitManager.h"




int main()
{
    Simulation simulation = Simulation(30);
    InitManager initManager = InitManager(&simulation);
    initManager.init_simulation(BreakingDam);
    // TestManager::test_kernel_integral();
    simulation.run();
    // simulation.render_from_file("testrender");
    return 0;
}
