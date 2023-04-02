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
	_stiffness = Parameters::GAMMA;
	_gravity.x = 0;
	_gravity.y = Parameters::GRAVITY;
	_viscosity = Parameters::VISCOSITY;

	_spawnLocations = std::vector<sf::Vector2f>();
	_spawnVelocities = std::vector<sf::Vector2f>();


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
	_averageDensity = 0;
	float density;
	float densityError;
	sf::Vector2f v_ij;
	sf::Vector2f x_ij;
	sf::Vector2f distance;
	sf::Vector2f viscosity = sf::Vector2f();
	sf::Vector2f c_f = sf::Vector2f();
	sf::Vector2f d_ij = sf::Vector2f();
	sf::Vector2f v_adv_ij;
	float distanceNorm = 1;
	float velocityNorm;
	sf::Vector2f kernelDeriv;
	float p_adv = 0;
	float p_0 = 0;
	float a_ii = 0;
	int l = 0;
	float Ap = 0;

	int numParticles = _particles.size();
	_numFluidParticles = 0;

	// Find Neighbors and count fluid particles
	for (int i = 0; i < numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		_numFluidParticles++;
		_particles[i]._neighbors = _hashManager.return_neighbors(&_particles[i], _neighborRadius);
	}
	if (_numFluidParticles == 0) { return; }

	// =========== PREDICT ADVECTION ================

	// Compute density, v_adv and c_i
	for (int i = 0; i < numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		density = 0;
		viscosity.x = 0;
		viscosity.y = 0;
		c_f.x = 0;
		c_f.y = 0;
		// Sum up density
		for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
			distanceNorm = Functions::calculate_distance_norm(Functions::calculate_distance(
				_particles[i]._position, _particles[i]._neighbors[j]->_position));
			if (_particles[i]._neighbors[j]->_type == fluid) {
				density += FluidParticle::_mass * Functions::kernel(distanceNorm);
			}
			else if(_particles[i]._neighbors[j]->_type == solid) {
				density += SolidParticle::_mass * Functions::kernel(distanceNorm);
			}
		}
		// Sum up viscosity and factor c
		// Viscosity need density of neighbors. Problem!
		for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
			x_ij = _particles[i]._position - _particles[i]._neighbors[j]->_position;
			v_ij = _particles[i]._velocity - _particles[i]._neighbors[j]->_velocity;
			distanceNorm = Functions::calculate_distance_norm(x_ij);
			kernelDeriv = Functions::kernel_derivation(x_ij, distanceNorm);
			if (_particles[i]._neighbors[j]->_type == fluid) {
				//viscosity += (FluidParticle::_mass * Functions::scalar_product2D(v_ij, x_ij))
				//	/ (_particles[i]._neighbors[j]->_density * Parameters::H * 0.01f
				//		+ Functions::scalar_product2D(x_ij, x_ij)
				//		* Parameters::H) * kernelDeriv;
				c_f += - (FluidParticle::_mass / (density * density)) * kernelDeriv;
			}
			else {
				// Uses rest density instead of boundary density atm. add factor gamma maybe?
				//viscosity += (SolidParticle::_mass * Functions::scalar_product2D(v_ij, x_ij))
				//	/ (FluidParticle::_restDensity * Parameters::H * 0.01f
				//		+ Functions::scalar_product2D(x_ij, x_ij)
				//		* Parameters::H) * Functions::kernel_derivation(distance, distanceNorm);
				c_f += - 2 * Parameters::GAMMA * (SolidParticle::_mass / (density * density)) * kernelDeriv;
			}
		}
		_particles[i]._density = density;
		if (_particles[i]._id == _watchedParticleId) { _watchedParticleDensity = density; }
		_averageDensity += std::max(density, FluidParticle::_restDensity);
		// Let's ignore viscosity for the time being
		//_particles[i]._v_adv = _particles[i]._velocity + Parameters::timeStepSize * _gravity
		//	+ Parameters::timeStepSize * FluidParticle::_materialParameter * 2 *
		//	(2 + distanceNorm) * viscosity;
		_particles[i]._v_adv = _particles[i]._velocity + Parameters::timeStepSize * _gravity;
		_particles[i].c_f = c_f;
	}
	// Division by _numFluidParticles not a problem, as function would have terminated if _nFP == 0
	_averageDensity = _averageDensity / _numFluidParticles;

	// Calculate predicted density, initial density and diagonal element a_ii
	for (int i = 0; i < numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		p_adv = 0;
		a_ii = 0;
		for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
			x_ij = _particles[i]._position - _particles[i]._neighbors[j]->_position;
			v_adv_ij = _particles[i]._v_adv - _particles[i]._neighbors[j]->_v_adv;
			distanceNorm = Functions::calculate_distance_norm(x_ij);
			kernelDeriv = Functions::kernel_derivation(x_ij, distanceNorm);

			if (_particles[i]._neighbors[j]->_type == fluid) {
				p_adv += FluidParticle::_mass
					* Functions::scalar_product2D(v_adv_ij, kernelDeriv);
				a_ii += FluidParticle::_mass *
					Functions::scalar_product2D(_particles[i].c_f, kernelDeriv);
				a_ii += + FluidParticle::_mass * Functions::scalar_product2D(
					(FluidParticle::_mass / (_particles[i]._density
						* _particles[i]._density)) * - kernelDeriv, kernelDeriv);
			}
			else if (_particles[i]._neighbors[j]->_type == solid) {
				p_adv += SolidParticle::_mass
					* Functions::scalar_product2D(v_adv_ij, kernelDeriv);
				a_ii += SolidParticle::_mass *
					Functions::scalar_product2D(_particles[i].c_f, kernelDeriv);
			}
		}
		_particles[i]._s_i = FluidParticle::_restDensity - _particles[i]._density - Parameters::timeStepSize * p_adv;
		_particles[i]._a_ii = Parameters::timeStepSize * Parameters::timeStepSize * a_ii;
		_particles[i]._pressure = std::max(Parameters::OMEGA * _particles[i]._s_i / _particles[i]._a_ii, 0.f);
		if (_particles[i]._a_ii == 0) { _particles[i]._pressure = 0; }
	}

	// ========================= Solve Pressure ============================
	while (true) {
		// Exit Condition
		if (l >= Parameters::MAX_SOLVER_ITERATIONS) {
			break;
		}

		for (int i = 0; i < numParticles; i++) {
			if (_particles[i]._type == solid) { continue; }
			_particles[i]._pressureAcc.x = 0;
			_particles[i]._pressureAcc.y = 0;
			for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
				x_ij = _particles[i]._position - _particles[i]._neighbors[j]->_position;
				distanceNorm = Functions::calculate_distance_norm(x_ij);
				kernelDeriv = Functions::kernel_derivation(x_ij, distanceNorm);

				if (_particles[i]._neighbors[j]->_type == fluid) {
					_particles[i]._pressureAcc += - FluidParticle::_mass * (_particles[i]._pressure
						/ (_particles[i]._density * _particles[i]._density)
						+ _particles[i]._neighbors[j]->_pressure
						/ (_particles[i]._neighbors[j]->_density
							* _particles[i]._neighbors[j]->_density))
						* kernelDeriv;
				}
				else {
					_particles[i]._pressureAcc += - Parameters::GAMMA * SolidParticle::_mass
						* (_particles[i]._pressure
						/ (_particles[i]._density * _particles[i]._density)) * 2
						* kernelDeriv;
				}
			}
		}
		densityError = 0;
		for (int i = 0; i < numParticles; i++) {
			if (_particles[i]._type == solid) { continue; }
			Ap = 0;
			for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
				x_ij = _particles[i]._position - _particles[i]._neighbors[j]->_position;
				distanceNorm = Functions::calculate_distance_norm(x_ij);
				kernelDeriv = Functions::kernel_derivation(x_ij, distanceNorm);
				if (_particles[i]._neighbors[j]->_type == fluid) {
					Ap += FluidParticle::_mass
						* Functions::scalar_product2D(_particles[i]._pressureAcc
							- _particles[i]._neighbors[j]->_pressureAcc, kernelDeriv);
				}
				else {
					Ap += SolidParticle::_mass
						* Functions::scalar_product2D(_particles[i]._pressureAcc, kernelDeriv);
				}
			}
			Ap = Ap * Parameters::timeStepSize * Parameters::timeStepSize;
			if (_particles[i]._a_ii != 0) {
				_particles[i]._pressure = std::max(_particles[i]._pressure + Parameters::OMEGA
					* (_particles[i]._s_i - Ap) / _particles[i]._a_ii, 0.f);
				if (Parameters::COLOR_CODE_PRESSURE) {
					_particles[i]._colorFactor = std::min((int)(_particles[i]._pressure * 0.005), 255);
				}
			}
			densityError += std::abs(Ap - _particles[i]._s_i);
			
		}
		// This is alright, _numFluidParticles is sure to be unequal to 0
		densityError = densityError / _numFluidParticles;

		// Increment Solver iteration
		l++;

		// Exit Condition
		if (densityError < Parameters::MAX_DENSITY_ERROR) {
			break;
		}
	}




	// =================Update Position and Velocity ========================
	if (!_moveParticles) { return; }
	for (int i = 0; i < numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		_particles[i]._velocity = Parameters::timeStepSize * _particles[i]._pressureAcc
			+ _particles[i]._v_adv;
		_particles[i]._position += Parameters::timeStepSize * _particles[i]._velocity;

		velocityNorm = Functions::calculate_distance_norm(_particles[i]._velocity);
		if (Parameters::COLOR_CODE_SPEED) {
			_particles[i]._colorFactor = std::min((int)(velocityNorm * 0.6), 255);
		}
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
void Simulation::spawn_particles() {
	if (_numIterations * Parameters::timeStepSize - _lastSpawnTime < _spawnDelay || _numFluidParticles >= maxNumParticles) { return; }
	for (int i = 0; i < _spawnLocations.size(); i++) {
		_particles.push_back(FluidParticle(_particles.size(), _spawnLocations[i], _spawnVelocities[i]));
	}
	_lastSpawnTime = _numIterations * Parameters::timeStepSize;
}



// _________________________________________________________________________________
void Simulation::run() {

	if (Parameters::INTERACTIVE) {
		_videoMode = sf::VideoMode(Parameters::WINDOW_WIDTH, Parameters::WINDOW_HEIGHT);
		_window.create(_videoMode, "SPH Fluid Solver");
		bool runSimulation = true;
		sf::Time newTime;
		sf::Time renderTime;
		float elapsedTime;
		_numUpdatesPerSec = 0;
		_numIterations = 0;
		_lastSpawnTime = 0.4;
		_renderer.init_solids(&_particles);
		_watchedParticleId = (int)(_particles.size() / 2);
		while (runSimulation) {
			// Update Clock
			newTime = _clock.getElapsedTime();
			elapsedTime = newTime.asSeconds() - _lastUpdate.asSeconds();
			_numUpdatesPerSec = 1 / elapsedTime;

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
				while (_particles[_watchedParticleId]._type != fluid) {
					int randomValue = std::rand() % _particles.size();
					_watchedParticleId = randomValue;
				}
			}

			if (_spawnParticles) {
				spawn_particles();
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
			_currentTime = _clock.getElapsedTime();

			_renderer.update_information(_currentTime.asSeconds(),
				_numIterations * Parameters::timeStepSize, _particles.size(),
				_numFluidParticles, _numUpdatesPerSec, _averageDensity, _maxVelocity,
				_watchedParticleDensity);

			_renderer.draw(&_window, &_particles, _watchedParticleId, _markedParticlesId, _testedParticlesId);

			renderTime += _clock.getElapsedTime() - _currentTime;

			_numIterations++;
		}
		std::cout << renderTime.asSeconds() / _numIterations << std::endl;

	}

	else {
		if (!Parameters::JUST_RENDER) {
			_numIterations = 0;
			_lastSpawnTime = 0;
			_renderFile.open("./renderFile.dat", std::fstream::out | std::fstream::trunc);
			_avgDensityFile.open("./avgDensityFile.dat", std::fstream::out | std::fstream::trunc);
			_renderFile << Parameters::timeStepSize << std::endl;
			for (int i = 0; i < _particles.size(); i++) {
				if (_particles[i]._type == solid) {
					_renderFile << "SP " << _particles[i]._position.x << " " << _particles[i]._position.y << " " <<
						0 << std::endl;
				}
			}

			_renderFile << "EOS" << std::endl;

			for (int timeSteps = 0; timeSteps < Parameters::SIMULATION_LENGTH; timeSteps++) {
				if (timeSteps % 10 == 0) {
					std::cout << "Calculating step " << timeSteps << " / " << Parameters::SIMULATION_LENGTH << "\r";
				}

				if (_spawnParticles) {
					spawn_particles();
				}
				update_hashTable();
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
				_numIterations++;
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
		int numSolids = 0;
		std::string type;
		float xValue;
		float yValue;
		float colorFactor;
		float density;
		bool runSimulation = true;
		float timeStepSize;
		int c = 0;
		_renderFile >> timeStepSize;
		_clock.restart();

		std::vector<sf::CircleShape> solidShapes = std::vector<sf::CircleShape>();

		while (true) {
			_renderFile >> type;
			if (type[0] == 'E') { break; }
			_renderFile >> xValue >> yValue >> colorFactor;
			sf::CircleShape newShape = sf::CircleShape();
			newShape.setFillColor(sf::Color::Color::White);
			newShape.setRadius(SolidParticle::_size * _zoomFactor);
			newShape.setPosition(xValue * _zoomFactor, yValue * _zoomFactor);
			solidShapes.push_back(newShape);
			numSolids++;
		}

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

				if (type[0] == 'F') {
					_renderer._fluidShape.setPosition(sf::Vector2f(xValue* _zoomFactor, yValue* _zoomFactor));
					_renderer._fluidShape.setFillColor(sf::Color::Blue + sf::Color::Color(0, colorFactor, 0));
					numFluids++;
					_window.draw(_renderer._fluidShape);
				}

				numShapes++;
				_renderFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			}

			for (int i = 0; i < numSolids; i++) {
				_window.draw(solidShapes[i]);
			}

 			_renderer.update_information(elapsedTime.asSeconds(), timeSteps* timeStepSize, numShapes, numFluids, numUpdatesPerSec, _averageDensity, _maxVelocity);
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