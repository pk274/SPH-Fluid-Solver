#pragma once
// Paul Kull, 2022
#include <vector>
#include <SFML/Graphics.hpp>

#include "./Particle.h"
#include "./SolidParticle.h"
#include "./FluidParticle.h"

#include "Renderer.h"

enum SimulationPreset {
	Empty = 0,
	StuffedBox = 1
};



class Simulation {
  public:
	  // std::vector<SolidParticle> _solidParticles;
	  // std::vector<FluidParticle> _fluidParticles;
	  std::vector<Particle> _particles;

	  Renderer renderer;
	  sf::VideoMode _videoMode;
	  sf::RenderWindow _window;

	  Simulation(SimulationPreset);
	  void init_empty_simulation();
	  void init_stuffed_box_simulation();

	  void update_physics();

	  void run();
};