// Paul Kull, 2022

#define _USE_MATH_DEFINES

#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "./Simulation.h"

constexpr double NEIGHBORHOOD_RADIUS = 5;
constexpr int NUM_SUPPOSED_NEIGHBORS = 7;
constexpr double GRAVITY = 9.8;
constexpr double STIFFNESS = 200;	// Has to go up with increasing mass of particles
constexpr double VISCOSITY = 2;	// Has to go down with increasing mass of particles
constexpr float H = 2.5;				// Distance of 2*H is supported by kernel -> H = neigRad / 2
constexpr float timeStepSize = 0.05;

// _________________________________________________________________________________
Simulation::Simulation(SimulationPreset preset = SmallBox, int framelimit) {

	_particles = std::vector<Particle>();
	_testedParticlesId = std::vector<int>();
	_particles.clear();

	_neighborRadius = NEIGHBORHOOD_RADIUS;
	_stiffness = STIFFNESS;
	_gravity.x = 0;
	_gravity.y = GRAVITY;
	_viscosity = VISCOSITY;


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
		init_breaking_dam_simulation(40, 10);
		break;
	case BigBreakingDam:
		init_breaking_dam_simulation(70, 5);
	}

	_clock = sf::Clock();
	_lastUpdate = _clock.getElapsedTime();
	std::srand(std::time(nullptr));
	_watchedParticleId = _particles.size() - 1;
	_hashManager = HashManager(_neighborRadius, 300);

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
void Simulation::init_empty_simulation() {
	
	_videoMode = sf::VideoMode(800, 800);
	_zoomFactor = 50;
	_window.create(_videoMode, "SPH Fluid Solver");
	_renderer = Renderer(_zoomFactor, FluidParticle::_size, SolidParticle::_size, _neighborRadius);

	sf::Vector2f pos = sf::Vector2f(5, 5);
	_particles.push_back(Particle(fluid, 0, pos));
	_moveParticles = false;
}


// _________________________________________________________________________________
void  Simulation::init_stuffed_box_simulation(int size, int zoom) {

	_videoMode = sf::VideoMode(800, 800);
	_window.create(_videoMode, "SPH Fluid Solver");
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
	_printFPS = false;
	_printParticleInfo = false;
}

// ______________________________________________________________________________________________________
void Simulation::init_single_particle_simulation(int size, int zoom) {

	_videoMode = sf::VideoMode(800, 800);
	_window.create(_videoMode, "SPH Fluid Solver");
	_zoomFactor = zoom;
	_renderer = Renderer(_zoomFactor, FluidParticle::_size, SolidParticle::_size, _neighborRadius);

	// Add Particles
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	pos.y = size * SolidParticle::_size;
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < size; j++) {
			_particles.push_back(SolidParticle(_particles.size(), pos));
			pos.x += SolidParticle::_size;
		}
		pos.x = 0;
		pos.y += SolidParticle::_size;
	}


	pos.x = 16;
	pos.y = 35;
	_particles.push_back(FluidParticle(_particles.size(), pos));

	_moveParticles = true;
	_testNeighbors = false;
	_testKernel = false;
	_printFPS = false;
	_printParticleInfo = true;
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
}

// _________________________________________________________________________________
void Simulation::init_random_particles_simulation(int size, int zoom, int numParticles) {
	
	_videoMode = sf::VideoMode(800, 800);
	_window.create(_videoMode, "SPH Fluid Solver");
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
	pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * 2);
	box = placeBox(pos, size - 4, _particles.size());
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size * 3, SolidParticle::_size * 3);
	box = placeBox(pos, size - 6, _particles.size());
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
}


// _________________________________________________________________________________
void  Simulation::init_breaking_dam_simulation(int size, int zoom) {
	
	_videoMode = sf::VideoMode(800, 800);
	_window.create(_videoMode, "SPH Fluid Solver");
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
	pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * 2);
	box = placeBox(pos, size - 4, _particles.size());
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size * 3, SolidParticle::_size * 3);
	box = placeBox(pos, size - 6, _particles.size());
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size * 4, SolidParticle::_size * (size - 4));
	for (int i = 0; i < size / 1.5; i++) {
		for (int j = 0; j < size / 3; j++) {
			_particles.push_back(FluidParticle(_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		pos.x = SolidParticle::_size * 4;
		pos.y -= FluidParticle::_size;
	}
	_moveParticles = true;
	_testNeighbors = false;
	_testKernel = false;
	_printFPS = false;
	_printParticleInfo = false;
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
	int numParticles = _particles.size();
	double density;
	float h = H;
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

	// Find Neighbors
	for (int i = 0; i < numParticles; i++) {
		_particles[i]._neighbors = _hashManager.return_neighbors(&_particles[i], _neighborRadius);
		if (_particles[i]._id == _watchedParticleId) {
			for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
				_markedParticlesId.push_back(_particles[i]._neighbors[j]->_id);
			}
		}
	}

	// Update Density and Pressure
	for (int i = 0; i < numParticles; i++) {
		density = 0;
		for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
			distanceNorm = Functions::calculate_distance_norm(Functions::calculate_distance(
				_particles[i]._position, _particles[i]._neighbors[j]->_position));
			if (_particles[i]._neighbors[j]->_type == fluid) {
				density += FluidParticle::_mass * Functions::kernel(distanceNorm, H);
			}
			else {
				float x = SolidParticle::_mass;
				float y = Functions::kernel(distanceNorm, H);
				density += SolidParticle::_mass * Functions::kernel(distanceNorm, H);
			}
		}
		_particles[i]._density = density;
		_particles[i]._pressure = std::max(0., _stiffness * (_particles[i]._density / FluidParticle::_restDensity - 1));
	}

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
					(Functions::scalar_product2D(x_ij, x_ij) + 0.01f * h * h)) *
					Functions::kernel_derivation(distance, distanceNorm, H);
				firstFraction = _particles[i]._pressure / (_particles[i]._density * _particles[i]._density);
				secondFraction = firstFraction;
				a_p -= (float)(firstFraction + secondFraction) *
					Functions::kernel_derivation(distance, distanceNorm, H) * SolidParticle::_mass;
			}
			else {
				a_nonp += (float)((FluidParticle::_mass * Functions::scalar_product2D(v_ij, x_ij)) /
					(Functions::scalar_product2D(x_ij, x_ij) + 0.01f * h * h)) *
					Functions::kernel_derivation(distance, distanceNorm, h);
				firstFraction = _particles[i]._pressure / (_particles[i]._density * _particles[i]._density);
				secondFraction = _particles[i]._neighbors[j]->_pressure /
					(_particles[i]._neighbors[j]->_density * _particles[i]._neighbors[j]->_density);
				a_p -= (float)(firstFraction + secondFraction) *
					Functions::kernel_derivation(distance, distanceNorm, h) * FluidParticle::_mass;
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
		_particles[i]._velocity += timeStepSize * _particles[i]._acceleration;
		_particles[i]._position += timeStepSize * _particles[i]._velocity;
		_particles[i]._lastUpdated = _clock.getElapsedTime();			// NOT necessary right now
	}
}




// _________________________________________________________________________________
void Simulation::run() {
	bool runSimulation = true;
	sf::Time elapsedTime;
	float numUpdatesPerSec = 0;

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
		}

		update_hashTable();
		update_physics();

		if (_testNeighbors) {
			_testedParticlesId.clear();
			_testedParticlesId = TestManager::test_correct_neighbor_amount(&_particles, NUM_SUPPOSED_NEIGHBORS);
		}
		if (_testKernel) {
			_testedParticlesId.clear();
			_testedParticlesId = TestManager::test_kernel(&_particles, H, _watchedParticleId);
		}
		if (_printFPS) {
			_renderer.update_information(_particles.size(), numUpdatesPerSec);
		}

		_renderer.update_graphics(&_particles, _watchedParticleId, _markedParticlesId, _testedParticlesId);
		_renderer.draw(&_window);
	}
}