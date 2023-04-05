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

#include "./Renderer.h"



class Simulation {
  public:
	  int _numIterations;

	  std::vector<Particle> _particles;
	  float _timeStepSize;
	  float _neighborRadius;
	  double _stiffness;
	  sf::Vector2f _gravity;
	  double _viscosity;
	  int _numParticles;
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
	  bool _spawnParticles;

	  float _spawnDelay;
	  float _lastSpawnTime;
	  int maxNumParticles;

	  bool _pauseSimulation;
	  bool _endSimulation;


	  sf::Clock _clock;
	  sf::Time _currentTime;
	  float _numUpdatesPerSec;
	  int _watchedParticleId;
	  std::vector<int> _markedParticlesId;
	  std::vector<int> _testedParticlesId;
	  sf::Time _lastUpdate;
	  float _averageDensity;
	  float _maxVelocity;
	  float _watchedParticleDensity;
	  float _avgNeighborhoodTime;
	  int _totalNumSolverIterations;
	  double _simulatedTime;
	  double _cflNumber;
	  double _nextFrame;
	  double _frameDistance;

	  std::vector<sf::Vector2f> _spawnLocations;
	  std::vector<sf::Vector2f> _spawnVelocities;

	  std::fstream _avgDensityFile;
	  std::fstream _renderFile;

	  Simulation(int framelimit = 30);

	  void check_input();
	  void run_tests();

	  void update_hashTable();
	  void update_hashTable_old();

	  void calculate_s_vd();
	  void calculate_s_di();

	  void jacobi_solve_vd();
	  void jacobi_solve_di();
	  void jacobi_solve_vd_ps();
	  void update_x_and_v();
	  void spawn_particles();
	  void delete_particles();

	  void update_physics();

	  void run();

	  void render_from_file(std::string fileName);
};