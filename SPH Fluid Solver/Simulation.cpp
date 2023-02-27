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
Simulation::Simulation(int framelimit) {

	_particles = std::vector<Particle>();
	_testedParticlesId = std::vector<int>();

	_neighborRadius = Parameters::NEIGHBORHOOD_RADIUS;
	_stiffness = Parameters::STIFFNESS;
	_gravity.x = 0;
	_gravity.y = Parameters::GRAVITY;
	_viscosity = Parameters::VISCOSITY;


	_clock = sf::Clock();
	_lastUpdate = _clock.getElapsedTime();
	std::srand(std::time(nullptr));
	_watchedParticleId = - 1;
	_maxVelocity = 0;

}



// _________________________________________________________________________________
void Simulation::update_hashTable() {
	_hashManager.update(_particles.size());
	for (int i = 0; i < _particles.size(); i++) {
		_hashManager.insert_item(&_particles[i]);
	}
}

// _________________________________________________________________________________
void Simulation::update_hashTable_old() {
	_hashManager.reset_buckets();
	for (int i = 0; i < _particles.size(); i++) {
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
	float velocityNorm;

	int numParticles = _particles.size();
	_numFluidParticles = 0;

	// Find Neighbors
	for (int i = 0; i < numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		// _particles[i]._neighbors = Functions::n_square_neighborhood_search(&_particles, i, _neighborRadius);
		_particles[i]._neighbors = _hashManager.return_neighbors(&_particles[i], _neighborRadius);
		if (_particles[i]._id == _watchedParticleId) {
			for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
				_markedParticlesId.push_back(_particles[i]._neighbors[j]->_id);
			}
		}
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
		if (_particles[i]._id == _watchedParticleId) { _watchedParticleDensity = _particles[i]._density; }
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

		velocityNorm = Functions::calculate_distance_norm(_particles[i]._velocity);
		_particles[i]._colorFactor = std::min((int)(velocityNorm * 0.6), 255);
		if (velocityNorm > _maxVelocity) {
			_maxVelocity = Functions::calculate_distance_norm(_particles[i]._velocity);
		}
	}

	// Delete stray particles
	if (!_deleteParticles) { return; }
	for (int i = 0; i < _particles.size(); i++) {
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
		sf::Time newTime;
		sf::Time renderTime;
		sf::Time currentTime;
		float elapsedTime;
		float numUpdatesPerSec = 0;
		int numIterations = 0;
		_renderer.init_solids(&_particles);
		while (runSimulation) {
			// Update Clock
			newTime = _clock.getElapsedTime();
			elapsedTime = newTime.asSeconds() - _lastUpdate.asSeconds();
			numUpdatesPerSec = 1 / elapsedTime;

			_lastUpdate = _clock.getElapsedTime();

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

			//update_hashTable_old();
			//update_hashTable();
			//update_physics();

			if (_testNeighbors) {
				_testedParticlesId.clear();
				_testedParticlesId = TestManager::test_correct_neighbor_amount(&_particles);
			}
			if (_testKernel) {
				_testedParticlesId.clear();
				_testedParticlesId = TestManager::test_kernel(&_particles, Parameters::H);
			}
			currentTime = _clock.getElapsedTime();

			_renderer.update_information(newTime.asSeconds(), _particles.size(),
				_numFluidParticles, numUpdatesPerSec, _averageDensity, _maxVelocity,
				_watchedParticleDensity);

			_renderer.draw(&_window, &_particles, _watchedParticleId, _markedParticlesId, _testedParticlesId);

			renderTime += _clock.getElapsedTime() - currentTime;

			numIterations++;
		}
		std::cout << renderTime.asSeconds() / numIterations << std::endl;

	}

	else {
		if (!Parameters::JUST_RENDER) {
			_renderFile.open("./renderFile.dat", std::fstream::out | std::fstream::trunc);
			_avgDensityFile.open("./avgDensityFile.dat", std::fstream::out | std::fstream::trunc);
			_renderFile << Parameters::timeStepSize << std::endl;
			for (int i = 0; i < _particles.size(); i++) {
				if (_particles[i]._type == solid) {
					_renderFile << "SP " << _particles[i]._position.x << " " << _particles[i]._position.y << " " <<
						0 << std::endl;
				}
			}

			for (int timeSteps = 0; timeSteps < Parameters::SIMULATION_LENGTH; timeSteps++) {
				if (timeSteps % 10 == 0) {
					std::cout << "Calculating step " << timeSteps << " / " << Parameters::SIMULATION_LENGTH << "\r";
				}
				update_hashTable();
				//update_hashTable_old();
				update_physics();
				int numParticles = _particles.size();
				if (timeSteps % Parameters::SPEEDUP == 0) {
					for (int i = 0; i < numParticles; i++) {
						if (_particles[i]._type == fluid) {
							_renderFile << "FP " << _particles[i]._position.x << " " << _particles[i]._position.y << " " <<
								_particles[i]._colorFactor << std::endl;
						}
					}
					_renderFile << "D " << _averageDensity << "\n";
					_renderFile << "EOU" << std::endl;
				}
	
				_avgDensityFile << timeSteps * Parameters::timeStepSize << " " << _averageDensity << "\n";
			}
			_renderFile.close();
		}
		_avgDensityFile.close();
		std::cout << _avgNeighborhoodTime / Parameters::SIMULATION_LENGTH << std::endl;

		_avgDensityFile.open("./avgDensityFile.dat", std::fstream::out | std::fstream::trunc);
		_renderFile.open("./renderFile.dat", std::fstream::in);
		_videoMode = sf::VideoMode(Parameters::WINDOW_WIDTH, Parameters::WINDOW_HEIGHT);
		_window.create(_videoMode, "SPH Fluid Solver");

		std::cout << "=======================================================" << std::endl;
		std::cout << "COMPUTATIONS READY! Press Space to watch the Simulation" << std::endl;
		std::cout << "=======================================================" << std::endl;

		while (!sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) {
			continue;
		}

		sf::Time elapsedTime;
		int numUpdatesPerSec;
		int numShapes = 0;
		int numFluids = 0;
		std::string type;
		float xValue;
		float yValue;
		float colorFactor;
		float density;
		bool runSimulation = true;
		float timeStepSize;
		int c = 0;
		_renderFile >> timeStepSize;

		for (int timeSteps = 0; timeSteps < Parameters::SIMULATION_LENGTH / Parameters::SPEEDUP; timeSteps++) {
			if (c % Parameters::RENDER_SPEEDUP != 0) {
				while (true) {
					_renderFile >> type;
					if (type[0] == 'E') { break; }
					_renderFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				}
				c++;
				continue;
			}
			elapsedTime = _clock.getElapsedTime();
			numShapes = 0;
			numFluids = 0;
			numUpdatesPerSec = 1 / (elapsedTime.asSeconds() - _lastUpdate.asSeconds());
			_lastUpdate = elapsedTime;
			_window.clear();
			

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
					_avgDensityFile << c * Parameters::timeStepSize << " " << _averageDensity << "\n";
					continue;
				}
				_renderFile >> xValue >> yValue >> colorFactor;
			
				
				if (type[0] == 'S') {
					_renderer._fluidShape.setPosition(sf::Vector2f(xValue* _zoomFactor, yValue* _zoomFactor));
					_renderer._fluidShape.setFillColor(sf::Color::White);
					_window.draw(_renderer._fluidShape);
				}
				else if (type[0] == 'F') {
					_renderer._fluidShape.setPosition(sf::Vector2f(xValue* _zoomFactor, yValue* _zoomFactor));
					_renderer._fluidShape.setFillColor(sf::Color::Blue + sf::Color::Color(0, colorFactor, 0));
					numFluids++;
					_window.draw(_renderer._fluidShape);
				}

				numShapes++;
				_renderFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			}

 			_renderer.update_information(timeSteps * timeStepSize, numShapes, numFluids, numUpdatesPerSec, _averageDensity, _maxVelocity);
			_window.draw(_renderer._infoPanel);
			_window.draw(_renderer._graphBackground);
			for (int i = 0; i < _renderer._graphShapes.size(); i++) {
				_window.draw(_renderer._graphShapes[i]);
			}
			_window.draw(_renderer._description);
			_window.draw(_renderer._information);
			_window.display();
			c++;
			if (!runSimulation) { break; }
		}
		_window.close();

	}

	_renderFile.close();
}



// ________________________________________________________________________
void Simulation::render_from_file(std::string fileName) {
	_renderFile.open("./simulationen/" + fileName + ".dat", std::fstream::in);

	std::cout << "=======================================================" << std::endl;
	std::cout << "COMPUTATIONS READY! Press Space to watch the Simulation" << std::endl;
	std::cout << "=======================================================" << std::endl;

	while (!sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) {
		continue;
	}
	_videoMode = sf::VideoMode(Parameters::WINDOW_WIDTH, Parameters::WINDOW_HEIGHT);
	_window.create(_videoMode, "SPH Fluid Solver");

	sf::Time elapsedTime;
	int numUpdatesPerSec;
	int numShapes = 0;
	int numFluids = 0;
	std::string type;
	float xValue;
	float yValue;
	float colorFactor;
	float density;
	bool runSimulation = true;
	float timeStepSize;
	int c = 0;
	_renderFile >> timeStepSize;

	for (int timeSteps = 0; timeSteps < Parameters::SIMULATION_LENGTH / Parameters::SPEEDUP; timeSteps++) {
		if (c % Parameters::RENDER_SPEEDUP != 0) {
			while (true) {
				_renderFile >> type;
				if (type[0] == 'E') { break; }
				_renderFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			}
			c++;
			continue;
		}
		elapsedTime = _clock.getElapsedTime();
		numShapes = 0;
		numFluids = 0;
		numUpdatesPerSec = 1 / (elapsedTime.asSeconds() - _lastUpdate.asSeconds());
		_lastUpdate = elapsedTime;
		_window.clear();


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
				_avgDensityFile << c * Parameters::timeStepSize << " " << _averageDensity << "\n";
				continue;
			}
			_renderFile >> xValue >> yValue >> colorFactor;


			if (type[0] == 'S') {
				_renderer._fluidShape.setPosition(sf::Vector2f(xValue * _zoomFactor, yValue * _zoomFactor));
				_renderer._fluidShape.setFillColor(sf::Color::White);
				_window.draw(_renderer._fluidShape);
			}
			else if (type[0] == 'F') {
				_renderer._fluidShape.setPosition(sf::Vector2f(xValue * _zoomFactor, yValue * _zoomFactor));
				_renderer._fluidShape.setFillColor(sf::Color::Blue + sf::Color::Color(0, colorFactor, 0));
				numFluids++;
				_window.draw(_renderer._fluidShape);
			}

			numShapes++;
			_renderFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		_renderer.update_information(timeSteps * timeStepSize, numShapes, numFluids, numUpdatesPerSec, _averageDensity, _maxVelocity);
		_window.draw(_renderer._infoPanel);
		_window.draw(_renderer._graphBackground);
		for (int i = 0; i < _renderer._graphShapes.size(); i++) {
			_window.draw(_renderer._graphShapes[i]);
		}
		_window.draw(_renderer._description);
		_window.draw(_renderer._information);
		_window.display();
		c++;
		if (!runSimulation) { break; }
	}
	_window.close();
	_renderFile.close();

}