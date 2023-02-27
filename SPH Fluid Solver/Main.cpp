// Paul Kull, 2022

#include <iostream>
#include <SFML/Graphics.hpp>

#include "./Simulation.h"
#include "./InitManager.h"

//____________________________________________
// TODO:
// -neighborhood search und hastable update verbessern
// _____________________________________________




int main()
{
    Simulation simulation = Simulation();
    InitManager initManager = InitManager(&simulation);
    initManager.init_simulation(GiantBox);
    // TestManager::test_kernel_integral();
    simulation.run();
    // simulation.render_from_file("testrender");
    return 0;
}
