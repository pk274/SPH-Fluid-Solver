#pragma once
// Paul Kull, 2022
#include <vector>
#include <SFML/Graphics.hpp>
#include <fstream>

#include "./Functions.h"
#include "./Particle.h"
#include "./SolidParticle.h"
#include "./FluidParticle.h"
#include "./HashManager.h"
#include "./TestManager.h"
#include "./CompactHashManager.h"

#include "Renderer.h"

enum SimulationPreset {
	Alone = 0,
	SmallBox = 1,
	StuffedBox = 2,
	SingleParticle = 3,
	RotatedBox = 4,
	FewParticles = 5,
	GiantBox = 6,
	FourLayers = 7,
	ManyLayers = 8,
	WideLayers = 9,
	BreakingDam = 10,
	BigBreakingDam = 11,
	Cup = 20,
	Complex = 30,
	Osmosis = 40,
};



class Simulation {
  public:
	  std::vector<Particle> _particles;
	  float _neighborRadius;
	  double _stiffness;
	  sf::Vector2f _gravity;
	  double _viscosity;
	  int _numFluidParticles;

	  Renderer _renderer;
	  sf::VideoMode _videoMode;
	  sf::RenderWindow _window;
	  int _zoomFactor;

	  CompactHashManager _hashManager;
	  // HashManager _hashManager;

	  bool _moveParticles;
	  bool _testNeighbors;
	  bool _testKernel;
	  bool _printFPS;
	  bool _printParticleInfo;
	  bool _deleteParticles;

	  sf::Clock _clock;
	  int _watchedParticleId;
	  std::vector<int> _markedParticlesId;
	  std::vector<int> _testedParticlesId;
	  sf::Time _lastUpdate;
	  float _averageDensity;
	  float _maxVelocity;
	  float _watchedParticleDensity;
	  float _avgNeighborhoodTime;

	  std::fstream _avgDensityFile;
	  std::fstream _renderFile;

	  Simulation(SimulationPreset preset, int framelimit = 30);
	  void init_empty_simulation();
	  void init_stuffed_box_simulation(int size = 50, int zoom = 5);
	  void init_single_particle_simulation(int size = 50, int zoom = 5);
	  void init_rotated_box_simulation(int size, int zoom, int rotation);
	  void init_random_particles_simulation(int size, int zoom, int numParticles);
	  void init_breaking_dam_simulation(int size = 50, int zoom = 5);
	  void init_layer_simulation(int size, int zoom, int layers);
	  void init_wide_layer_simulation(int size, int zoom, int layers);
	  void init_cup_simulation(int size, int zoom);
	  void init_complex_simulation(int size, int zoom);
	  void init_osmosis_simulation(int size, int zoom);

	  void update_hashTable();
	  void update_hashTable_old();

	  void update_physics();

	  void run();

	  void render_from_file(std::string fileName);
};