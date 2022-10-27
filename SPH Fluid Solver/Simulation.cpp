// Paul Kull, 2022

#define _USE_MATH_DEFINES

#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>

#include "./Simulation.h"

// _________________________________________________________________________________
Simulation::Simulation(SimulationPreset preset = StuffedBox, int framelimit) {

	_particles = std::vector<Particle>();
	_particles.clear();

	_neighborRadius = 5;
	_stiffness = 1;
	_gravity.x = 0;
	_gravity.y = -9.8;
	_viscosity = 1;

	switch (preset) {
	case Empty:
		init_empty_simulation();
		break;
	case SmallBox:
		init_small_box();
		break;
	case StuffedBox:
		init_stuffed_box_simulation(35, 11);
		break;
	case StuffedBoxZoomed:
		init_stuffed_box_zoomed_simulation();
		break;
	}

	_clock = sf::Clock();
	_lastUpdate = _clock.getElapsedTime();
	std::srand(std::time(nullptr));
	_watchedParticleId = -100;
	_hashManager = HashManager(_neighborRadius);

}

// _________________________________________________________________________________
std::vector<Particle> placeBox(sf::Vector2f pos, int size) {
	std::vector<Particle> box = std::vector<Particle>();
	for (int i = 0; i < size; i++) {
		box.push_back(SolidParticle(box.size(), pos));
		pos.y += 2;
	}
	for (int i = 0; i < size; i++) {
		box.push_back(SolidParticle(box.size(), pos));
		pos.x += 2;
	}
	for (int i = 0; i < size; i++) {
		box.push_back(SolidParticle(box.size(), pos));
		pos.y -= 2;
	}
	for (int i = 0; i < size; i++) {
		box.push_back(SolidParticle(box.size(), pos));
		pos.x -= 2;
	}
	return box;
}


// _________________________________________________________________________________
void Simulation::init_empty_simulation() {
	sf::VideoMode _videoMode = sf::VideoMode();
	_videoMode.size = sf::Vector2u(200, 200);
	_window.create(_videoMode, "SPH Fluid Solver");
	_renderer = Renderer(1);
}

// _________________________________________________________________________________
void Simulation::init_small_box() {
	sf::VideoMode _videoMode = sf::VideoMode();
	_videoMode.size = sf::Vector2u(800, 800);
	_window.create(_videoMode, "SPH Fluid Solver");
	_zoomFactor = 3;
	_renderer = Renderer(_zoomFactor, FluidParticle::_size, SolidParticle::_size, _neighborRadius);


	sf::Vector2f pos = sf::Vector2f(150, 200);
	std::vector<Particle> box = placeBox(pos, 4);
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}

	pos.x = 152;
	pos.y = 202;

	for (int i = 0; i < 3; i++) {
		for (int ii = 0; ii < 3; ii++) {
			_particles.push_back(FluidParticle(_particles.size(), pos));
			pos.y += 2;
		}
		pos.y = 202;
		pos.x += 2;
	}

}


// _________________________________________________________________________________
void  Simulation::init_stuffed_box_simulation(int size, int zoom) {

	sf::VideoMode _videoMode = sf::VideoMode();
	_videoMode.size = sf::Vector2u(800, 800);
	_window.create(_videoMode, "SPH Fluid Solver");
	_zoomFactor = zoom;
	_renderer = Renderer(_zoomFactor, FluidParticle::_size, SolidParticle::_size, _neighborRadius);

	// Add Particles
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_particles.push_back(box[i]);
	}


	pos.x = 2;
	pos.y = 2;

	for (int i = 0; i < size - 1; i++) {
		for (int ii = 0; ii < size - 1; ii++) {
			_particles.push_back(FluidParticle(_particles.size(), pos));
			pos.y += 2;
		}
		pos.y = 2;
		pos.x += 2;
	}
}


// _________________________________________________________________________________
void  Simulation::init_stuffed_box_zoomed_simulation() {

	sf::VideoMode _videoMode = sf::VideoMode();
	_videoMode.size = sf::Vector2u(800, 800);
	_window.create(_videoMode, "SPH Fluid Solver");
	_zoomFactor = 10;
	_renderer = Renderer(_zoomFactor, FluidParticle::_size, SolidParticle::_size, _neighborRadius);

	// Add Particles
	sf::Vector2f pos = sf::Vector2f(150, 200);
	for (int i = 0; i < 25; i++) {
		_particles.push_back(SolidParticle(_particles.size(), pos));
		pos.y += 20;
	}
	for (int i = 0; i < 25; i++) {
		_particles.push_back(SolidParticle(_particles.size(), pos));
		pos.x += 20;
	}
	for (int i = 0; i < 25; i++) {
		_particles.push_back(SolidParticle(_particles.size(), pos));
		pos.y -= 20;
	}
	for (int i = 0; i < 25; i++) {
		_particles.push_back(SolidParticle(_particles.size(), pos));
		pos.x -= 20;
	}

	pos.x = 170;
	pos.y = 220;

	for (int i = 0; i < 24; i++) {
		for (int ii = 0; ii < 24; ii++) {
			_particles.push_back(FluidParticle(_particles.size(), pos));
			pos.y += 20;
		}
		pos.y = 220;
		pos.x += 20;
	}
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
float kernel(float distance, int h = 1) {
	float q = distance / h;
	float a = 5 / (14 * M_PI * pow(h, 2));
	if (distance < 1) {
		return a * (pow(2 - q, 3) - 4 * pow(1-q, 3));
	}
	if (distance < 2) {
		return a * (pow(2 - q, 3));
	}
	return 0;
}

// _________________________________________________________________________________
sf::Vector2f kernel_derivation(sf::Vector2f distance, float distanceNorm, int h = 1) {
	float q = distanceNorm / h;
	if (distanceNorm < 1) {
		return sf::Vector2f((float)(15 * q + (3 * q - 4) / 14 * pow(h, 2) * M_PI) * distance);
	}
	if (distanceNorm < 2) {
		return sf::Vector2f((float)(- 15 * pow(q - 2, 2) / 14 * pow(h, 2) * M_PI) * distance);
	}
	return sf::Vector2f(0, 0);
}




// _________________________________________________________________________________
// ONLY WORKS IF THERE ARE ONLY FLUID PARTICLES, BECAUSE OF FLUIDPARTICLE::_MASS
void Simulation::update_physics() {
	_markedParticlesId.clear();
	int numParticles = _particles.size();
	std::vector<Particle*> neighbors = std::vector<Particle*>();
	float density;
	int h = 1;
	sf::Vector2f a_nonp;
	sf::Vector2f a_p;
	sf::Vector2f a;
	sf::Vector2f v_ij;
	for (int i = 0; i < numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		neighbors = _hashManager.return_neighbors(&_particles[i], _neighborRadius);
		if (_particles[i]._id == _watchedParticleId) {
			for (int j = 0; j < neighbors.size(); j++) {
				_markedParticlesId.push_back(neighbors[j]->_id);
			}
		}
		density = 0;
		for (int j = 0; j < neighbors.size(); j++) {
			density += FluidParticle::_mass * kernel(neighbors[j]->_distanceNorm);
		}
		_particles[i]._density = density;
		_particles[i]._pressure = std::max(0.f, _stiffness * (_particles[i]._density / FluidParticle::_restDensity - 1));
		a_nonp.x = 0;
		a_nonp.y = 0;
		for (int j = 0; j < neighbors.size(); j++) {
			v_ij = _particles[i]._velocity - neighbors[j]->_velocity;
			a_nonp += (FluidParticle::_mass * v_ij * neighbors[j]->_distanceNorm) /
				(float)(pow(neighbors[j]->_distanceNorm, 2) + 0.01 * h * h) * kernel_derivation(neighbors[j]->_distanceNorm);
		}
		a_nonp *= 2 * _viscosity;
		a_nonp += _gravity;
		a_p.x = 0;
		a_p.y = 0;
		for (int j = 0; j < neighbors.size(); j++) {
			a_p -= _particles[i]._pressure / pow(_particles[i]._density, 2) + neighbors[j]->_pressure / pow(neighbors[j]->_density, 2)
					* kernel_derivation(neighbors[j]->_distance) * FluidParticle::_mass;
		}
		
	}
}




// _________________________________________________________________________________
void Simulation::run() {
	bool runSimulation = true;
	sf::Time elapsedTime;
	float numUpdatesPerSec = 0;

	while (runSimulation) {
		elapsedTime = _clock.getElapsedTime();
		numUpdatesPerSec = 1 / (elapsedTime.asSeconds() - _lastUpdate.asSeconds());
		_lastUpdate = elapsedTime;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
			_window.close();
			runSimulation = false;
			break;
		}
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::R)) {
			int randomValue = std::rand() % _particles.size();
			_watchedParticleId = randomValue;
		}
		update_hashTable();
		update_physics();
		_renderer.update_graphics(&_particles, _watchedParticleId, _markedParticlesId);
		_renderer.update_information(_particles.size(), numUpdatesPerSec);
		_renderer.draw(&_window);
	}
}