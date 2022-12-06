// Paul Kull, 2022

#define _USE_MATH_DEFINES

#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>

#include "./Simulation.h"
#include "Parameters.h"


// _________________________________________________________________________________
Simulation::Simulation(SimulationPreset preset = SmallBox, int framelimit) {

	_particles = std::vector<Particle>();
	_testedParticlesId = std::vector<int>();
	_particles.clear();

	_neighborRadius = Parameters::NEIGHBORHOOD_RADIUS;
	_stiffness = Parameters::STIFFNESS;
	_gravity.x = 0;
	_gravity.y = Parameters::GRAVITY;
	_viscosity = Parameters::VISCOSITY;


	switch (preset) {
	case Alone:
		init_empty_simulation();
		break;
	case SmallBox:
		init_stuffed_box_simulation(20, 17);
		break;
	case StuffedBox:
		init_stuffed_box_simulation(35, 11);
		break;
	case SingleParticle:
		init_single_particle_simulation(35, 10);
		break;
	case RotatedBox:
		init_rotated_box_simulation(20, 8, 30);
		break;
	case FewParticles:
		init_random_particles_simulation(35, 10, 10);
		break;
	case BreakingDam:
		init_breaking_dam_simulation(39, 10);
		break;
	case BigBreakingDam:
		init_breaking_dam_simulation(76, 5);
		break;
	case GiantFuckingBox:
		init_stuffed_box_simulation(300, 1);
		break;
	case FourLayers:
		init_layer_simulation(39, 10, 4);
		break;
	case ManyLayers:
		init_layer_simulation(39, 10, 10);
		break;
	}

	_clock = sf::Clock();
	_lastUpdate = _clock.getElapsedTime();
	std::srand(std::time(nullptr));
	_watchedParticleId = _particles.size() - 1;
	_hashManager = HashManager(_neighborRadius, 300);
	_maxTimeStep = 0;

	_avgDensityFile.open("./avgDensityFile", std::fstream::out | std::fstream::trunc);

}

// _________________________________________________________________________________
std::vector<Particle> placeBox(sf::Vector2f pos, int size, int startId = 0) {
	std::vector<Particle> box = std::vector<Particle>();
	for (int i = 0; i < size; i++) {
		box.push_back(SolidParticle(startId + box.size(), pos));
		pos.y += SolidParticle::_size;
	}
	for (int i = 0; i < size; i++) {
		box.push_back(SolidParticle(startId + box.size(), pos));
		pos.x += SolidParticle::_size;
	}
	for (int i = 0; i < size; i++) {
		box.push_back(SolidParticle(startId + box.size(), pos));
		pos.y -= SolidParticle::_size;
	}
	for (int i = 0; i < size; i++) {
		box.push_back(SolidParticle(startId + box.size(), pos));
		pos.x -= SolidParticle::_size;
	}
	return box;
}

// _________________________________________________________________________________
std::vector<Particle> placeTallBox(sf::Vector2f pos, int size, int minus = 0, int startId = 0) {
	std::vector<Particle> box = std::vector<Particle>();
	for (int i = 0; i < size; i++) {
		box.push_back(SolidParticle(startId + box.size(), pos));
		pos.y += SolidParticle::_size;
	}
	for (int i = 0; i < size / 2 - minus; i++) {
		box.push_back(SolidParticle(startId + box.size(), pos));
		pos.x += SolidParticle::_size;
	}
	for (int i = 0; i < size; i++) {
		box.push_back(SolidParticle(startId + box.size(), pos));
		pos.y -= SolidParticle::_size;
	}
	for (int i = 0; i < size / 2 - minus; i++) {
		box.push_back(SolidParticle(startId + box.size(), pos));
		pos.x -= SolidParticle::_size;
	}
	return box;
}


// _________________________________________________________________________________
void Simulation::init_empty_simulation() {
	
	_videoMode = sf::VideoMode(Parameters::WINDOW_WIDTH, Parameters::WINDOW_HEIGHT);
	_zoomFactor = 50;
	_window.create(_videoMode, "SPH Fluid Solver");
	_renderer = Renderer(_zoomFactor, FluidParticle::_size, SolidParticle::_size, _neighborRadius);

	sf::Vector2f pos = sf::Vector2f(5, 5);
	_particles.push_back(Particle(fluid, 0, pos));
	_moveParticles = false;
}


// _________________________________________________________________________________
void  Simulation::init_stuffed_box_simulation(int size, int zoom) {

	_zoomFactor = zoom;
	_renderer = Renderer(_zoomFactor, FluidParticle::_size, SolidParticle::_size, _neighborRadius);

	// Add Particles
	// Add Particles for arena
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeBox(pos, size - 2, _particles.size());
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * 2);
	box = placeBox(pos, size - 4, _particles.size());
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}

	pos = sf::Vector2f(SolidParticle::_size * 3, SolidParticle::_size * 3);

	for (int i = 0; i < size - 5; i++) {
		for (int ii = 0; ii < size - 5; ii++) {
			_particles.push_back(FluidParticle(_particles.size(), pos));
			pos.y += FluidParticle::_size;
		}
		pos.y = SolidParticle::_size * 3;
		pos.x += FluidParticle::_size;
	}
	_moveParticles = false;
	_testNeighbors = false;
	_testKernel = true;
	_printFPS = true;
	_printParticleInfo = false;
	_deleteParticles = false;
}

// ______________________________________________________________________________________________________
void Simulation::init_single_particle_simulation(int size, int zoom) {

	_zoomFactor = zoom;
	_renderer = Renderer(_zoomFactor, FluidParticle::_size, SolidParticle::_size, _neighborRadius);

	// Add Particles
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeBox(pos, size - 2, _particles.size());
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}


	pos.x = 16;
	pos.y = 35;
	_particles.push_back(FluidParticle(_particles.size(), pos));

	_moveParticles = true;
	_testNeighbors = false;
	_testKernel = false;
	_printFPS = true;
	_printParticleInfo = true;
	_deleteParticles = true;
}

// _________________________________________________________________________________
void Simulation::init_rotated_box_simulation(int size, int zoom, int rotation) {
	
	init_stuffed_box_simulation(size, zoom);
	sf::Vector2f offset = sf::Vector2f(0, size * 1.5);
	for (int i = 0; i < _particles.size(); i++) {
		_particles[i]._position.x += std::cos(rotation * 3.141592653589 / 180) * _particles[i]._position.x
			+ std::sin(rotation * 3.141592653589 / 180) * _particles[i]._position.y;
		_particles[i]._position.y += - std::sin(rotation * 3.141592653589 / 180) * _particles[i]._position.x
			+ std::cos(rotation * 3.141592653589 / 180) * _particles[i]._position.y;
		_particles[i]._position += offset;
	}
	_testNeighbors = true;
	_testKernel = false;
	_printFPS = true;
	_printParticleInfo = false;
	_deleteParticles = false;
}

// _________________________________________________________________________________
void Simulation::init_random_particles_simulation(int size, int zoom, int numParticles) {
	
	_zoomFactor = zoom;
	_renderer = Renderer(_zoomFactor, FluidParticle::_size, SolidParticle::_size, _neighborRadius);

	// Add Particles for arena
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeBox(pos, size - 2, _particles.size());
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	for (int i = 0; i < numParticles; i++) {
		float x = (std::rand() % (size - 10)) + SolidParticle::_size * 5;
		float y = (std::rand() % (size - 10)) + SolidParticle::_size * 5;
		_particles.push_back(FluidParticle(_particles.size(), sf::Vector2f(x, y)));
	}
	_moveParticles = true;
	_testNeighbors = false;
	_testKernel = false;
	_printFPS = false;
	_printParticleInfo = false;
	_deleteParticles = true;
}


// _________________________________________________________________________________
void  Simulation::init_breaking_dam_simulation(int size, int zoom) {
	
	_zoomFactor = zoom;
	_renderer = Renderer(_zoomFactor, FluidParticle::_size, SolidParticle::_size, _neighborRadius);

	// Add Particles for arena
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeBox(pos, size - 2, _particles.size());
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * (size - 2));
	for (int i = 0; i < size / 1.5; i++) {
		for (int j = 0; j < size / 3; j++) {
			_particles.push_back(FluidParticle(_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		pos.x = SolidParticle::_size * 2;
		pos.y -= FluidParticle::_size;
	}
	_moveParticles = true;
	_testNeighbors = false;
	_testKernel = false;
	_printFPS = true;
	_printParticleInfo = false;
	_deleteParticles = true;
}

// _________________________________________________________________________________
void Simulation::init_layer_simulation(int size, int zoom, int layers) {
	_zoomFactor = zoom;
	_renderer = Renderer(_zoomFactor, FluidParticle::_size, SolidParticle::_size, _neighborRadius);

	// Add Particles for arena
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeTallBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeTallBox(pos, size - 2, 1, _particles.size());
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * (size - 2));
	for (int i = 0; i < layers; i++) {
		for (int j = 0; j < size / 2 - 3; j++) {
			_particles.push_back(FluidParticle(_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		pos.x = SolidParticle::_size * 2;
		pos.y -= FluidParticle::_size;
	}
	_moveParticles = true;
	_testNeighbors = false;
	_testKernel = false;
	_printFPS = true;
	_printParticleInfo = false;
	_deleteParticles = false;
}



// _________________________________________________________________________________
void Simulation::update_hashTable() {
	_hashManager.reset_buckets();
	int numParticles = _particles.size();
	for (int i = 0; i < numParticles; i++) {
		_hashManager.insert_item(&_particles[i]);
	}
}


// _________________________________________________________________________________
void Simulation::update_physics() {
	_markedParticlesId.clear();
	_averageDensity = 0;
	double density;
	float h = Parameters::H;
	float particleMass;
	sf::Vector2f a_nonp;
	sf::Vector2f a_p;
	sf::Vector2f a;
	sf::Vector2f v_ij;
	sf::Vector2f x_ij;
	sf::Vector2f distance;
	float distanceNorm;
	double firstFraction;
	double secondFraction;

	int numParticles = _particles.size();
	_numFluidParticles = 0;

	// Find Neighbors
	for (int i = 0; i < numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		_particles[i]._neighbors = _hashManager.return_neighbors(&_particles[i], _neighborRadius);
		// if (_particles[i]._id == _watchedParticleId) {
		// 	for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
		// 		_markedParticlesId.push_back(_particles[i]._neighbors[j]->_id);
		// 	}
		// }
	}

	// Update Density and Pressure
	for (int i = 0; i < numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		_numFluidParticles++;
		density = 0;
		for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
			distanceNorm = Functions::calculate_distance_norm(Functions::calculate_distance(
				_particles[i]._position, _particles[i]._neighbors[j]->_position));
			if (_particles[i]._neighbors[j]->_type == fluid) {
				density += FluidParticle::_mass * Functions::kernel(distanceNorm);
			}
			else {
				density += SolidParticle::_mass * Functions::kernel(distanceNorm);
			}
		}
		_particles[i]._density = density;
		_averageDensity += std::max((float) density, FluidParticle::_restDensity);
		_particles[i]._pressure = std::max(0., _stiffness * (_particles[i]._density / FluidParticle::_restDensity - 1));
	}
	_averageDensity = _averageDensity / _numFluidParticles;

	// Update Accelerations
	for (int i = 0; i < numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		a_nonp.x = 0;
		a_nonp.y = 0;
		a_p.x = 0;
		a_p.y = 0;
		for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
			if (_particles[i]._neighbors[j]->_id == _particles[i]._id) { continue; }
			distance = Functions::calculate_distance(_particles[i]._position, _particles[i]._neighbors[j]->_position);
			distanceNorm = Functions::calculate_distance_norm(distance);
			v_ij = _particles[i]._velocity - _particles[i]._neighbors[j]->_velocity;
			x_ij = _particles[i]._position - _particles[i]._neighbors[j]->_position;

			if (_particles[i]._neighbors[j]->_type == solid) {
				a_nonp += (float)((SolidParticle::_mass * Functions::scalar_product2D(v_ij, x_ij)) /
					(Functions::scalar_product2D(x_ij, x_ij) + 0.01f * h * h))
					* Functions::kernel_derivation(distance, distanceNorm);
				firstFraction = _particles[i]._pressure / (_particles[i]._density * _particles[i]._density);
				secondFraction = firstFraction;
				a_p -= (float)(firstFraction + secondFraction) *
					Functions::kernel_derivation(distance, distanceNorm) * SolidParticle::_mass;
			}
			else {
				a_nonp += (float)((FluidParticle::_mass * Functions::scalar_product2D(v_ij, x_ij)) /
					(Functions::scalar_product2D(x_ij, x_ij) + 0.01f * h * h
					* _particles[i]._neighbors[j]->_density))
					* Functions::kernel_derivation(distance, distanceNorm);
				firstFraction = _particles[i]._pressure / (_particles[i]._density * _particles[i]._density);
				secondFraction = _particles[i]._neighbors[j]->_pressure /
					(_particles[i]._neighbors[j]->_density * _particles[i]._neighbors[j]->_density);
				a_p -= (float)(firstFraction + secondFraction) *
					Functions::kernel_derivation(distance, distanceNorm) * FluidParticle::_mass;
			}
		}
		a_nonp *= (float)(2 * _viscosity);
		a_nonp += _gravity;
		_particles[i]._pressureAcc = a_p;
		_particles[i]._acceleration = a_p + a_nonp;
	}

	// Update Position and Velocity
	if (!_moveParticles) { return; }
	for (int i = 0; i < numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		_particles[i]._velocity += Parameters::timeStepSize * _particles[i]._acceleration;
		_particles[i]._position += Parameters::timeStepSize * _particles[i]._velocity;

		if (Functions::calculate_distance_norm(_particles[i]._velocity) > _maxVelocity) {
			_maxVelocity = Functions::calculate_distance_norm(_particles[i]._velocity);
			_maxTimeStep = 0.1 * Parameters::H / _maxVelocity;
		}
	}

	// Delete stray particles
	if (!_deleteParticles) { return; }
	for (int i = 0; i < _particles.size(); i++) {
		if (_particles[i]._id == _watchedParticleId) { _watchedParticleDensity = _particles[i]._density; }
		if (_particles[i]._type == solid) { continue; }
		// Delete Particles which fell out of the window
		if (_particles[i]._position.x < -10 || _particles[i]._position.y < -10 ||
			_particles[i]._position.x > Parameters::WINDOW_WIDTH / _zoomFactor ||
			_particles[i]._position.y > Parameters::WINDOW_HEIGHT / _zoomFactor) {
			_particles.erase(_particles.begin() + i);
			_numFluidParticles--;
		}
	}
}



// _________________________________________________________________________________
void Simulation::run() {

	if (Parameters::INTERACTIVE) {
		_videoMode = sf::VideoMode(Parameters::WINDOW_WIDTH, Parameters::WINDOW_HEIGHT);
		_window.create(_videoMode, "SPH Fluid Solver");
		bool runSimulation = true;
		sf::Time elapsedTime;
		float numUpdatesPerSec = 0;
		int numIterations = 0;
		while (runSimulation) {
			// Update Clock
			elapsedTime = _clock.getElapsedTime();
			numUpdatesPerSec = 1 / (elapsedTime.asSeconds() - _lastUpdate.asSeconds());

			_lastUpdate = elapsedTime;

			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) {
				// Pause
          		continue;
			}

			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
				// Close application
				_window.close();
				runSimulation = false;
				break;
			}
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::R)) {
				// Switch observed particle randomly
				int randomValue = std::rand() % _particles.size();
				_watchedParticleId = randomValue;
				//while (_particles[_watchedParticleId]._type != fluid) {
				//	int randomValue = std::rand() % _particles.size();
				//	_watchedParticleId = randomValue;
				//}
			}

			update_hashTable();
			update_physics();

			if (_testNeighbors) {
				_testedParticlesId.clear();
				_testedParticlesId = TestManager::test_correct_neighbor_amount(&_particles);
			}
			if (_testKernel) {
				_testedParticlesId.clear();
				_testedParticlesId = TestManager::test_kernel(&_particles, Parameters::H);
			}
			_renderer.update_information(elapsedTime.asSeconds(), _particles.size(),
				_numFluidParticles, numUpdatesPerSec, _averageDensity, _maxTimeStep,
				_watchedParticleDensity);

			_renderer.update_graphics(&_particles, _numFluidParticles, _watchedParticleId, _markedParticlesId, _testedParticlesId);
			_renderer.draw(&_window);

			_avgDensityFile << numIterations * Parameters::timeStepSize << " " << _averageDensity << "\n";
			numIterations++;
		}

	}

	else {
		if (!Parameters::JUST_RENDER) {
			_renderFile.open("./renderFile", std::fstream::out | std::fstream::trunc);
			_renderFile << Parameters::timeStepSize << std::endl;
			for (int i = 0; i < _particles.size(); i++) {
				if (_particles[i]._type == solid) {
					_renderFile << "SolidParticle " << _particles[i]._position.x << " " << _particles[i]._position.y << " " <<
						0 << std::endl;
				}
			}

			for (int timeSteps = 0; timeSteps < Parameters::SIMULATION_LENGTH; timeSteps++) {
				if (timeSteps % 10 == 0) {
					std::cout << "Calculating step " << timeSteps << " / " << Parameters::SIMULATION_LENGTH << "\r";
				}
				update_hashTable();
				update_physics();
				int numParticles = _particles.size();
				if (timeSteps % Parameters::SPEEDUP == 0) {
					for (int i = 0; i < numParticles; i++) {
						if (_particles[i]._type == fluid) {
							_renderFile << "FLuidParticle " << _particles[i]._position.x << " " << _particles[i]._position.y << " " <<
								std::min((int)(Functions::calculate_distance_norm(_particles[i]._velocity) * 0.6), 255) << std::endl;
						}
					}
					_renderFile << "Density " << _averageDensity << "\n";
					_renderFile << "END_OF_UPDATE" << std::endl;
				}
				_avgDensityFile << timeSteps * Parameters::timeStepSize << " " << _averageDensity << "\n";
			}
			_renderFile.close();
		}
		_avgDensityFile.close();

		std::cout << "=======================================================" << std::endl;
		std::cout << "COMPUTATIONS READY! Press Space to watch the Simulation" << std::endl;
		std::cout << "=======================================================" << std::endl;

		while (!sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) {
			continue;
		}

		_renderFile.open("./renderFile", std::fstream::in);
		// _avgDensityFile.open("./avgDensityFile", std::fstream::out | std::fstream::trunc);
		_videoMode = sf::VideoMode(Parameters::WINDOW_WIDTH, Parameters::WINDOW_HEIGHT);
		_window.create(_videoMode, "SPH Fluid Solver");
		sf::Time elapsedTime;
		int numUpdatesPerSec;
		int numShapes = 0;
		int numFluids = 0;
		std::string type;
		float xValue;
		float yValue;
		float speed;
		float density;
		bool runSimulation = true;
		float timeStepSize;
		_renderFile >> timeStepSize;

		for (int timeSteps = 0; timeSteps < Parameters::SIMULATION_LENGTH / Parameters::SPEEDUP; timeSteps++) {
			elapsedTime = _clock.getElapsedTime();
			numShapes = 0;
			numFluids = 0;
			numUpdatesPerSec = 1 / (elapsedTime.asSeconds() - _lastUpdate.asSeconds());
			_lastUpdate = elapsedTime;

			while (true) {
				if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
					// Close application
					runSimulation = false;
					break;
				}

				_renderFile >> type;
				if (type[0] == 'E') {
					break;
				}
				if (type[0] == 'D') {
					_renderFile >> density;
					_averageDensity = density;
					// _avgDensityFile << timeSteps * Parameters::timeStepSize << " " << _averageDensity << "\n";
					continue;
				}
				_renderFile >> xValue >> yValue >> speed;
			
				if (_renderer._fluidParticleShapes.size() + _renderer._solidParticleShapes.size() <= numShapes) {
					if (type[0] == 'S') {
						_renderer._solidParticleShapes.push_back(sf::CircleShape(_renderer._solidShapeRadius));
						_renderer._solidParticleShapes.back().setPosition(sf::Vector2f(xValue* _zoomFactor, yValue* _zoomFactor));
						_renderer._solidParticleShapes.back().setFillColor(sf::Color::White);
					}
					else if (type[0] == 'F') {
						_renderer._fluidParticleShapes.push_back(sf::CircleShape(_renderer._fluidShapeRadius));
						_renderer._fluidParticleShapes.back().setPosition(sf::Vector2f(xValue* _zoomFactor, yValue* _zoomFactor));
						_renderer._fluidParticleShapes.back().setFillColor(sf::Color::Blue + sf::Color::Color(0, speed, 0));
						numFluids++;
					}
				}
				else {
					if (type[0] == 'F') {
						_renderer._fluidParticleShapes.at(numFluids).setPosition(sf::Vector2f(xValue* _zoomFactor, yValue* _zoomFactor));
						_renderer._fluidParticleShapes.at(numFluids).setFillColor(FluidParticle::_stasisColor + sf::Color::Color(0, speed, 0));
						numFluids++;
					}
				}
				numShapes++;
				_renderFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			}
			while (numFluids < _renderer._fluidParticleShapes.size()) {
				_renderer._fluidParticleShapes.pop_back();
			}
			numShapes = numFluids + _renderer._solidParticleShapes.size();
 			_renderer.update_information(timeSteps * timeStepSize, numShapes, numFluids, numUpdatesPerSec, _averageDensity);
			_renderer.draw(&_window);
			if (!runSimulation) { break; }
		}
		_window.close();

	}

	// _avgDensityFile.close();
	_renderFile.close();
}



// ________________________________________________________________________
void Simulation::render_from_file(std::string fileName) {
	_renderFile.open("./renderFile", std::fstream::in);
	int timeStepSize;
	_renderFile >> timeStepSize;
	std::string type;
	float x;
	float y;
	float color;
	int numSolids = 0;
	int numFluids = 0;

	sf::Time elapsedTime;
	float numUpdatesPerSec = 0;

	// Add solid particles
	while (true) {
		_renderFile >> type;
		if (type[0] != 'S') { break; }
		_renderFile >> x >> y;
		_particles.push_back(SolidParticle(numSolids, sf::Vector2f(x, y)));
		numSolids++;
	}

	for (int i = 0; i < Parameters::SIMULATION_LENGTH; i++) {
		elapsedTime = _clock.getElapsedTime();
		numUpdatesPerSec = 1 / (elapsedTime.asSeconds() - _lastUpdate.asSeconds());
		while (true) {
			numFluids = 0;
			_renderFile >> type;
			if (type[0] == 'F') {
				_renderFile >> x >> y >> color;
				if (_particles.size() >= numSolids + numFluids) {
					_particles.at(numSolids + numFluids)._position = sf::Vector2f(x, y);
					numFluids++;
					continue;
				}
				else {
					_particles.push_back(FluidParticle(numFluids + numSolids, sf::Vector2f(x, y)));
					numFluids++;
					continue;
				}
			}
			if (type[0] == 'D') {
				_renderFile >> _averageDensity;
				continue;
			}
			if (type[0] == 'E') {
				break;
			}
		}
		while (numSolids + numFluids < _particles.size()) {
			_particles.pop_back();
		}
		_renderer.update_graphics(&_particles,numFluids, -1, _markedParticlesId, _testedParticlesId);
		_renderer.update_information(i * timeStepSize, numSolids + numFluids, numFluids, numUpdatesPerSec, _averageDensity, -1, -1);
		_renderer.draw(&_window);
	}
}