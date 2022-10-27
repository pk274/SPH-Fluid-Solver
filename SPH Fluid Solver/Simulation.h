#pragma once
// Paul Kull, 2022
#include <vector>
#include <SFML/Graphics.hpp>

#include "./Particle.h"
#include "./SolidParticle.h"
#include "./FluidParticle.h"
#include "./HashManager.h"

#include "Renderer.h"

enum SimulationPreset {
	Empty = 0,
	SmallBox = 1,
	StuffedBox = 2,
	StuffedBoxZoomed = 3
};



class Simulation {
  public:
	  // std::vector<SolidParticle> _solidParticles;
	  // std::vector<FluidParticle> _fluidParticles;
	  std::vector<Particle> _particles;
	  int _neighborRadius;
	  float _stiffness;
	  sf::Vector2f _gravity;
	  float _viscosity;

	  Renderer _renderer;
	  sf::VideoMode _videoMode;
	  sf::RenderWindow _window;
	  int _zoomFactor;

	  HashManager _hashManager;

	  sf::Clock _clock;
	  int _watchedParticleId;
	  std::vector<int> _markedParticlesId;
	  sf::Time _lastUpdate;

	  Simulation(SimulationPreset preset, int framelimit = 30);
	  void init_empty_simulation();
	  void init_small_box();
	  void init_stuffed_box_simulation(int size = 50, int zoom = 5);
	  void init_stuffed_box_zoomed_simulation();

	  void update_hashTable();
	  std::vector<Particle*> find_neighbors(Particle* particle);

	  void update_physics();

	  void run();
};