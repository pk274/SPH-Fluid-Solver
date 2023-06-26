// Paul Kull, 2022

#define _USE_MATH_DEFINES

#include <cstdlib>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>

#include "./Simulation.h"
#include "Parameters.h"


// _________________________________________________________________________________
Simulation::Simulation(int framesPerSec) {

	_particles = std::vector<Particle>();
	_testedParticlesId = std::vector<int>();

	_timeStepSize = Parameters::TIME_STEP;
	_simulatedTime = 0;

	_neighborRadius = Parameters::NEIGHBORHOOD_RADIUS;
	_stiffness = Parameters::GAMMA;
	_gravity.x = 0;
	_gravity.y = Parameters::GRAVITY;
	_viscosity = Parameters::BOUNDARY_VISCOSITY;

	_spawnLocations = std::vector<sf::Vector2f>();
	_spawnVelocities = std::vector<sf::Vector2f>();

	_pauseSimulation = false;
	_endSimulation = false;

	_clock = sf::Clock();
	_lastUpdate = _clock.getElapsedTime();
	std::srand(std::time(nullptr));
	_watchedParticleId = - 1;
	_maxVelocity = 1;
	_nextFrame = 0;
	_frameDistance = (int)10000 * (1. / framesPerSec);
	_frameDistance = _frameDistance / 10000;
	_moveSolids = false;

}

// _________________________________________________________________________________
void Simulation::check_input() {
	if (sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) {
		_pauseSimulation = true;
	}

	if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
		// Close application
		_window.close();
		_endSimulation = true;
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
void Simulation::EOS_solve() {
	_stiffness = Parameters::EOS_STIFFNESS;
	_viscosity = Parameters::EOS_VISCOSITY;
	float density;
	sf::Vector2f v_ij;
	sf::Vector2f x_ij;
	sf::Vector2f distance;
	sf::Vector2f a_p = sf::Vector2f();
	sf::Vector2f a_nonp = sf::Vector2f();
	float distanceNorm = 1;
	sf::Vector2f kernelDeriv;
	float kernel;
	float firstFraction = 0;
	float secondFraction = 0;
	float velocityNorm = 0;
	// Update Density and Pressure
	for (int i = 0; i < _numParticles; i++) {
		if (_particles[i]._type != fluid) { continue; }
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
		_averageDensity += std::max((float)density, FluidParticle::_restDensity);
		_particles[i]._pressure = std::max(0., _stiffness * (std::pow(_particles[i]._density / FluidParticle::_restDensity, 7) - 1));
	}
	_averageDensity = _averageDensity / _numFluidParticles;

	// Update Accelerations
	for (int i = 0; i < _numParticles; i++) {
		if (_particles[i]._id == _watchedParticleId) { _watchedParticleDensity = _particles[i]._density; }
		if (_particles[i]._type != fluid) { continue; }
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
					(Functions::scalar_product2D(x_ij, x_ij) + 0.01f * Parameters::H * Parameters::H))
					* Functions::kernel_derivation(distance, distanceNorm);
				firstFraction = _particles[i]._pressure / (_particles[i]._density * _particles[i]._density);
				secondFraction = firstFraction;
				a_p -= (float)(firstFraction + secondFraction) *
					Functions::kernel_derivation(distance, distanceNorm) * SolidParticle::_mass;
			}
			else {
				a_nonp += (float)((FluidParticle::_mass * Functions::scalar_product2D(v_ij, x_ij)) /
					(Functions::scalar_product2D(x_ij, x_ij) + 0.01f * Parameters::H * Parameters::H
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
		_particles[i]._v_adv = a_nonp;
		_particles[i]._pressureAcc = a_p;
	}

	// Update Position and Velocity
	if (!_moveParticles) { return; }
	_maxVelocity = 0;
	for (int i = 0; i < _numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		_particles[i]._velocity += _timeStepSize * (_particles[i]._pressureAcc + _particles[i]._v_adv);
		_particles[i]._position += _timeStepSize * _particles[i]._velocity;

		velocityNorm = Functions::calculate_distance_norm(_particles[i]._velocity);
		_particles[i]._colorFactor = std::min((int)(velocityNorm * 0.6), 255);
		if (velocityNorm > _maxVelocity) {
			_maxVelocity = Functions::calculate_distance_norm(_particles[i]._velocity);
		}
	}
}

// _________________________________________________________________________________
void Simulation::calculate_s_di() {
	float density;
	sf::Vector2f v_ij;
	sf::Vector2f x_ij;
	sf::Vector2f distance;
	sf::Vector2f viscosity = sf::Vector2f();
	sf::Vector2f c_f = sf::Vector2f();
	sf::Vector2f d_ij = sf::Vector2f();
	sf::Vector2f v_adv_ij;
	float distanceNorm = 1;
	sf::Vector2f kernelDeriv;
	float p_adv = 0;
	float p_0 = 0;
	float a_ii = 0;


	// =========== PREDICT ADVECTION ================

	// Compute density, v_adv and c_i
	for (int i = 0; i < _numParticles; i++) {
		if (_particles[i]._type != fluid) { continue; }
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
			else if (_particles[i]._neighbors[j]->_type == solid || _particles[i]._neighbors[j]->_type == moving) {
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
				viscosity += (FluidParticle::_mass * Functions::scalar_product2D(v_ij, x_ij))
					/ (_particles[i]._neighbors[j]->_density * (Parameters::H
						* Parameters::H * 0.01f
						+ Functions::scalar_product2D(x_ij, x_ij))) * kernelDeriv;
				c_f += -(FluidParticle::_mass / (density * density)) * kernelDeriv;
			}
			else {
				// Uses rest density instead of boundary density atm. add factor gamma maybe?
				viscosity += (SolidParticle::_mass * Functions::scalar_product2D(v_ij, x_ij))
					/ (density * (Parameters::H * 0.01f * Parameters::H
						+ Functions::scalar_product2D(x_ij, x_ij))) * Parameters::BOUNDARY_VISCOSITY * kernelDeriv;
				c_f += -2 * Parameters::GAMMA * (SolidParticle::_mass / (density * density)) * kernelDeriv;
			}
		}
		_particles[i]._density = density;
		if (_particles[i]._id == _watchedParticleId) { _watchedParticleDensity = density; }
		_averageDensity += std::max(density, FluidParticle::_restDensity);
		_particles[i]._v_adv = _particles[i]._velocity + _timeStepSize * (_gravity
			+ FluidParticle::_materialParameter * 8 * viscosity);
		_particles[i].c_f = c_f;
	}
	// Division by _numFluidParticles not a problem, as function would have terminated if _nFP == 0
	_averageDensity = _averageDensity / _numFluidParticles;

	// Calculate source term and diagonal element a_ii
	for (int i = 0; i < _numParticles; i++) {
		if (_particles[i]._type != fluid) { continue; }
		p_adv = 0;
		a_ii = 0;
		for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
			x_ij = _particles[i]._position - _particles[i]._neighbors[j]->_position;
			distanceNorm = Functions::calculate_distance_norm(x_ij);
			kernelDeriv = Functions::kernel_derivation(x_ij, distanceNorm);

			if (_particles[i]._neighbors[j]->_type == fluid) {
				v_adv_ij = _particles[i]._v_adv - _particles[i]._neighbors[j]->_v_adv;
				p_adv += FluidParticle::_mass
					* Functions::scalar_product2D(v_adv_ij, kernelDeriv);
				a_ii += FluidParticle::_mass *
					Functions::scalar_product2D(_particles[i].c_f, kernelDeriv);
				a_ii += FluidParticle::_mass * Functions::scalar_product2D(
					(FluidParticle::_mass / (_particles[i]._density
						* _particles[i]._density)) * -kernelDeriv, kernelDeriv);
			}
			else if (_particles[i]._neighbors[j]->_type == solid) {
				p_adv += SolidParticle::_mass
					* Functions::scalar_product2D(_particles[i]._v_adv, kernelDeriv);
				a_ii += SolidParticle::_mass *
					Functions::scalar_product2D(_particles[i].c_f, kernelDeriv);
			}
			else if (_particles[i]._neighbors[j]->_type == moving) {
				v_adv_ij = _particles[i]._v_adv - _particles[i]._neighbors[j]->_velocity;
				p_adv += SolidParticle::_mass
					* Functions::scalar_product2D(_particles[i]._v_adv, kernelDeriv);
				a_ii += SolidParticle::_mass *
					Functions::scalar_product2D(_particles[i].c_f, kernelDeriv);
			}
		}
		_particles[i]._s_i = FluidParticle::_restDensity - _particles[i]._density - _timeStepSize * p_adv;
		_particles[i]._a_ii = _timeStepSize * _timeStepSize * a_ii;
		_particles[i]._pressure = std::max(Parameters::OMEGA * _particles[i]._s_i / _particles[i]._a_ii, 0.f);
		if (_particles[i]._a_ii == 0) { _particles[i]._pressure = 0; }
	}
}

// _________________________________________________________________________________
void Simulation::calculate_s_vd() {
	float density;
	sf::Vector2f v_ij;
	sf::Vector2f x_ij;
	sf::Vector2f distance;
	sf::Vector2f viscosity = sf::Vector2f();
	sf::Vector2f c_f = sf::Vector2f();
	sf::Vector2f d_ij = sf::Vector2f();
	sf::Vector2f v_adv_ij;
	float distanceNorm = 1;
	sf::Vector2f kernelDeriv;
	float velocityDiv = 0;
	float p_0 = 0;
	float a_ii = 0;


	// =========== PREDICT ADVECTION ================

	// Compute density, v_adv and c_i
	for (int i = 0; i < _numParticles; i++) {
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
			else if (_particles[i]._neighbors[j]->_type == solid) {
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
				viscosity += (FluidParticle::_mass * Functions::scalar_product2D(v_ij, x_ij))
					/ (_particles[i]._neighbors[j]->_density * (Parameters::H
						* Parameters::H * 0.01f
						+ Functions::scalar_product2D(x_ij, x_ij))) * kernelDeriv;
				c_f += -(FluidParticle::_mass / (density * density)) * kernelDeriv;
			}
			else {
				// Uses rest density instead of boundary density atm. add factor gamma maybe?
				viscosity += (SolidParticle::_mass * Functions::scalar_product2D(v_ij, x_ij))
					/ (density * (Parameters::H * 0.01f * Parameters::H
						+ Functions::scalar_product2D(x_ij, x_ij))) * Parameters::BOUNDARY_VISCOSITY * kernelDeriv;
				c_f += -2 * Parameters::GAMMA * (SolidParticle::_mass / (density * density)) * kernelDeriv;
			}
		}
		_particles[i]._density = density;
		if (_particles[i]._id == _watchedParticleId) { _watchedParticleDensity = density; }
		_averageDensity += std::max(density, FluidParticle::_restDensity);
		_particles[i]._v_adv = _particles[i]._velocity + _timeStepSize * (_gravity
			+ FluidParticle::_materialParameter * 8 * viscosity);
		_particles[i].c_f = c_f;
	}
	// Division by _numFluidParticles not a problem, as function would have terminated if _nFP == 0
	_averageDensity = _averageDensity / _numFluidParticles;

	// Calculate source term and diagonal element a_ii
	for (int i = 0; i < _numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		velocityDiv = 0;
		a_ii = 0;
		for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
			x_ij = _particles[i]._position - _particles[i]._neighbors[j]->_position;
			distanceNorm = Functions::calculate_distance_norm(x_ij);
			kernelDeriv = Functions::kernel_derivation(x_ij, distanceNorm);

			if (_particles[i]._neighbors[j]->_type == fluid) {
				v_adv_ij = _particles[i]._v_adv - _particles[i]._neighbors[j]->_v_adv;
				velocityDiv -= (FluidParticle::_mass / _particles[i]._neighbors[j]->_density)
					* Functions::scalar_product2D(v_adv_ij, kernelDeriv);
				a_ii += FluidParticle::_mass *
					Functions::scalar_product2D(_particles[i].c_f, kernelDeriv);
				a_ii += FluidParticle::_mass * Functions::scalar_product2D(
					(FluidParticle::_mass / (_particles[i]._density
						* _particles[i]._density)) * -kernelDeriv, kernelDeriv);
			}
			else if (_particles[i]._neighbors[j]->_type == solid) {
				velocityDiv -= (SolidParticle::_mass / FluidParticle::_restDensity)
					* Functions::scalar_product2D(_particles[i]._v_adv, kernelDeriv);
				a_ii += SolidParticle::_mass *
					Functions::scalar_product2D(_particles[i].c_f, kernelDeriv);
			}
		}
		_particles[i]._s_i = _timeStepSize * _particles[i]._density * velocityDiv;
		_particles[i]._a_ii = _timeStepSize * _timeStepSize * a_ii;
		_particles[i]._pressure = std::max(Parameters::OMEGA * _particles[i]._s_i / _particles[i]._a_ii, 0.f);
		if (_particles[i]._a_ii == 0) { _particles[i]._pressure = 0; }
	}
}



// _________________________________________________________________________________
void Simulation::jacobi_solve() {

	int l = 0;
	float Ap = 0;
	sf::Vector2f x_ij;
	float distanceNorm;
	sf::Vector2f kernelDeriv;

	_numSolverIterations = 0;
	while (true) {
		// Exit Condition
		if (l >= Parameters::MAX_SOLVER_ITERATIONS) {
			break;
		}
		for (int i = 0; i < _numParticles; i++) {
			if (_particles[i]._type != fluid) { continue; }
			_particles[i]._pressureAcc.x = 0;
			_particles[i]._pressureAcc.y = 0;
			for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
				x_ij = _particles[i]._position - _particles[i]._neighbors[j]->_position;
				distanceNorm = Functions::calculate_distance_norm(x_ij);
				kernelDeriv = Functions::kernel_derivation(x_ij, distanceNorm);

				if (_particles[i]._neighbors[j]->_type == fluid) {
					_particles[i]._pressureAcc += -FluidParticle::_mass * (_particles[i]._pressure
						/ (_particles[i]._density * _particles[i]._density)
						+ _particles[i]._neighbors[j]->_pressure
						/ (_particles[i]._neighbors[j]->_density
							* _particles[i]._neighbors[j]->_density))
						* kernelDeriv;
				}
				else {
					_particles[i]._pressureAcc += - Parameters::GAMMA * SolidParticle::_mass
						* 2 * (_particles[i]._pressure
							/ (_particles[i]._density * _particles[i]._density))
						* kernelDeriv;
				}
			}
		}
		_estimatedDensityError = 0;
		int numParticlesWithNeighbors = 0;
		for (int i = 0; i < _numParticles; i++) {
			if (_particles[i]._type != fluid) { continue; }
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
			Ap = Ap * _timeStepSize * _timeStepSize;
			if (_particles[i]._a_ii != 0) {
				_particles[i]._pressure = std::max(_particles[i]._pressure + Parameters::OMEGA
					* (_particles[i]._s_i - Ap) / _particles[i]._a_ii, 0.f);
				//_estimatedDensityError += std::abs(Ap - _particles[i]._s_i);
				//numParticlesWithNeighbors++;
			}
			_estimatedDensityError += std::max(Ap - _particles[i]._s_i, 0.f);
			//_estimatedDensityError += Ap - _particles[i]._s_i;

		}
		// This is alright, _numFluidParticles is sure to be unequal to 0
		_estimatedDensityError = _estimatedDensityError / _numFluidParticles;
		//_estimatedDensityError = _estimatedDensityError / numParticlesWithNeighbors;

		// Increment Solver iteration
		l++;

		// Exit Condition
		if (_estimatedDensityError < Parameters::MAX_DENSITY_ERROR && l >= Parameters::MIN_SOLVER_ITERATIONS) {
			break;
		}
	}
	_numSolverIterations = l;
	_totalNumSolverIterations += _numSolverIterations;
	
}


// _________________________________________________________________________________
void Simulation::update_x_and_v() {
	float velocityNorm;
	_maxVelocity = 0;
	for (int i = 0; i < _numParticles; i++) {
		if (_particles[i]._type == solid) { continue; }
		if (_particles[i]._type == moving) {
			_particles[i]._position += _timeStepSize * _particles[i]._velocity;
			continue;
		}
		_particles[i]._velocity = _timeStepSize * _particles[i]._pressureAcc
			+ _particles[i]._v_adv;
		_particles[i]._position += _timeStepSize * _particles[i]._velocity;

		velocityNorm = Functions::calculate_distance_norm(_particles[i]._velocity);
		if (Parameters::COLOR_CODE_SPEED) {
			_particles[i]._colorFactor = std::min((int)(velocityNorm / Parameters::VELOCITY_CODE_RANGE), 255);
		}
		if (velocityNorm > _maxVelocity) {
			_maxVelocity = Functions::calculate_distance_norm(_particles[i]._velocity);
		}
	}
}

// _________________________________________________________________________________
void Simulation::update_moving_objects() {
	bool switchState;
	for (int i = 0; i < _movingObjects.size(); i++) {
		switchState = false;
		for (int ii = 0; ii < _movingObjects[i]._particles.size(); ii++) {
			if (_movingObjects[i]._conditionBigger[_movingObjects[i]._state]) {
				if (_movingObjects[i]._particles[ii]->_position.x >
					_movingObjects[i]._conditions[_movingObjects[i]._state].x
					||
					_movingObjects[i]._particles[ii]->_position.y >
					_movingObjects[i]._conditions[_movingObjects[i]._state].y) {
					switchState = true;
					break;
				}
			}
			else if (_movingObjects[i]._particles[ii]->_position.x <
				_movingObjects[i]._conditions[_movingObjects[i]._state].x
				||
				_movingObjects[i]._particles[ii]->_position.y <
				_movingObjects[i]._conditions[_movingObjects[i]._state].y) {
				switchState = true;
				break;
			}
		}
		if (switchState) {
			_movingObjects[i]._state++;
			if (_movingObjects[i]._state >= _movingObjects[i]._directions.size()) {
				_movingObjects[i]._state = 0;
			}
			for (int ii = 0; ii < _movingObjects[i]._particles.size(); ii++) {
				_movingObjects[i]._particles[ii]->_velocity
					= _movingObjects[i]._directions[_movingObjects[i]._state];
			}
		}
	}
}
// _________________________________________________________________________________
void Simulation::spawn_particles() {
	if (_simulatedTime - _lastSpawnTime < _spawnDelay || _numFluidParticles >= _maxNumParticles) { return; }
	for (int i = 0; i < _spawnLocations.size(); i++) {
		_particles.push_back(FluidParticle(_particles.size(), _spawnLocations[i], _spawnVelocities[i]));
	}
	_lastSpawnTime = _simulatedTime;
}

// _________________________________________________________________________________
void Simulation::delete_particles() {
	for (int i = 0; i < _particles.size(); i++) {
		if (_particles[i]._type != fluid) { continue; }
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
void Simulation::update_physics() {
	_averageDensity = 0;
	
	if (Parameters::SOLVE_PPE) {
		// Calculate source term and solve pressure
	    if (Parameters::S_VD) { calculate_s_vd(); jacobi_solve(); }
		else if (Parameters::S_DI){ calculate_s_di(); jacobi_solve(); }
		else if (Parameters::S_VD_DI) {
			calculate_s_vd();
			jacobi_solve();
			_averageDensity = 0;
			for (int i = 0; i < _numParticles; i++) {
				if (_particles[i]._type != fluid) { continue; }
				_particles[i]._velocity = _timeStepSize * _particles[i]._pressureAcc
					+ _particles[i]._v_adv;
				_particles[i]._position += _timeStepSize * _particles[i]._velocity;
			}
			float density;
			sf::Vector2f v_ij;
			sf::Vector2f x_ij;
			sf::Vector2f distance;
			sf::Vector2f c_f = sf::Vector2f();
			sf::Vector2f d_ij = sf::Vector2f();
			sf::Vector2f v_adv_ij;
			float distanceNorm = 1;
			sf::Vector2f kernelDeriv;
			float p_adv = 0;
			float p_0 = 0;
			float a_ii = 0;

			// Compute density, v_adv and c_i
			for (int i = 0; i < _numParticles; i++) {
				if (_particles[i]._type != fluid) { continue; }
				density = 0;
				c_f.x = 0;
				c_f.y = 0;
				// Sum up density
				for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
					distanceNorm = Functions::calculate_distance_norm(Functions::calculate_distance(
						_particles[i]._position, _particles[i]._neighbors[j]->_position));
					if (_particles[i]._neighbors[j]->_type == fluid) {
						density += FluidParticle::_mass * Functions::kernel(distanceNorm);
					}
					else if (_particles[i]._neighbors[j]->_type == solid || _particles[i]._neighbors[j]->_type == moving) {
						density += SolidParticle::_mass * Functions::kernel(distanceNorm);
					}
				}
				// Sum up c
				for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
					x_ij = _particles[i]._position - _particles[i]._neighbors[j]->_position;
					v_ij = _particles[i]._velocity - _particles[i]._neighbors[j]->_velocity;
					distanceNorm = Functions::calculate_distance_norm(x_ij);
					kernelDeriv = Functions::kernel_derivation(x_ij, distanceNorm);
					if (_particles[i]._neighbors[j]->_type == fluid) {
						c_f += -(FluidParticle::_mass / (density * density)) * kernelDeriv;
					}
					else {
						c_f += -2 * Parameters::GAMMA * (SolidParticle::_mass / (density * density)) * kernelDeriv;
					}
				}
				_particles[i]._density = density;
				if (_particles[i]._id == _watchedParticleId) { _watchedParticleDensity = density; }
				_averageDensity += std::max(density, FluidParticle::_restDensity);
				_particles[i]._v_adv = _particles[i]._velocity;
				_particles[i].c_f = c_f;
			}
			// Division by _numFluidParticles not a problem, as function would have terminated if _nFP == 0
			_averageDensity = _averageDensity / _numFluidParticles;

			// Calculate source term and diagonal element a_ii
			for (int i = 0; i < _numParticles; i++) {
				if (_particles[i]._type != fluid) { continue; }
				p_adv = 0;
				a_ii = 0;
				for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
					x_ij = _particles[i]._position - _particles[i]._neighbors[j]->_position;
					distanceNorm = Functions::calculate_distance_norm(x_ij);
					kernelDeriv = Functions::kernel_derivation(x_ij, distanceNorm);

					if (_particles[i]._neighbors[j]->_type == fluid) {
						v_adv_ij = _particles[i]._v_adv - _particles[i]._neighbors[j]->_v_adv;
						p_adv += FluidParticle::_mass
							* Functions::scalar_product2D(v_adv_ij, kernelDeriv);
						a_ii += FluidParticle::_mass *
							Functions::scalar_product2D(_particles[i].c_f, kernelDeriv);
						a_ii += FluidParticle::_mass * Functions::scalar_product2D(
							(FluidParticle::_mass / (_particles[i]._density
								* _particles[i]._density)) * -kernelDeriv, kernelDeriv);
					}
					else if (_particles[i]._neighbors[j]->_type == solid) {
						p_adv += SolidParticle::_mass
							* Functions::scalar_product2D(_particles[i]._v_adv, kernelDeriv);
						a_ii += SolidParticle::_mass *
							Functions::scalar_product2D(_particles[i].c_f, kernelDeriv);
					}
					else if (_particles[i]._neighbors[j]->_type == moving) {
						v_adv_ij = _particles[i]._v_adv - _particles[i]._neighbors[j]->_velocity;
						p_adv += SolidParticle::_mass
							* Functions::scalar_product2D(_particles[i]._v_adv, kernelDeriv);
						a_ii += SolidParticle::_mass *
							Functions::scalar_product2D(_particles[i].c_f, kernelDeriv);
					}
				}
				_particles[i]._s_i = FluidParticle::_restDensity - _particles[i]._density - _timeStepSize * p_adv;
				_particles[i]._a_ii = _timeStepSize * _timeStepSize * a_ii;
				_particles[i]._pressure = std::max(Parameters::OMEGA * _particles[i]._s_i / _particles[i]._a_ii, 0.f);
				if (_particles[i]._a_ii == 0) { _particles[i]._pressure = 0; }
			}
			jacobi_solve();
			//sf::Vector2f v_grad1;
			//sf::Vector2f v_grad2;
			//sf::Vector2f v_ji;
			//for (int i = 0; i < _numParticles; i++) {
			//	if (_particles[i]._type != fluid) { continue; }
			//	v_grad1.x = 0;
			//	v_grad1.y = 0;
			//	v_grad2.x = 0;
			//	v_grad2.y = 0;
			//	for (int j = 0; j < _particles[i]._neighbors.size(); j++) {
			//		x_ij = _particles[i]._position - _particles[i]._neighbors[j]->_position;
			//		distanceNorm = Functions::calculate_distance_norm(x_ij);
			//		kernelDeriv = Functions::kernel_derivation(x_ij, distanceNorm);
			//		if (_particles[i]._neighbors[j]->_type == fluid) {
			//			v_ji = _particles[i]._neighbors[j]->_velocity - _particles[i]._velocity;
			//			v_grad1.x += FluidParticle::_mass * (v_ji.x * kernelDeriv.x);
			//			v_grad1.y += FluidParticle::_mass * (v_ji.y * kernelDeriv.x);
			//			v_grad2.x += FluidParticle::_mass * (v_ji.x * kernelDeriv.y);
			//			v_grad2.y += FluidParticle::_mass * (v_ji.y * kernelDeriv.y);
			//		}
			//		else if (_particles[i]._neighbors[j]->_type == solid) {
			//			v_grad1.x += SolidParticle::_mass * (-_particles[i]._velocity.x * kernelDeriv.x);
			//			v_grad1.y += SolidParticle::_mass * (-_particles[i]._velocity.y * kernelDeriv.x);
			//			v_grad2.x += SolidParticle::_mass * (-_particles[i]._velocity.x * kernelDeriv.y);
			//			v_grad2.y += SolidParticle::_mass * (-_particles[i]._velocity.y * kernelDeriv.y);
			//		}
			//		else if (_particles[i]._neighbors[j]->_type == moving) {
			//			v_ji = _particles[i]._neighbors[j]->_velocity - _particles[i]._velocity;
			//			v_grad1.x += SolidParticle::_mass * (v_ji.x * kernelDeriv.x);
			//			v_grad1.y += SolidParticle::_mass * (v_ji.y * kernelDeriv.x);
			//			v_grad2.x += SolidParticle::_mass * (v_ji.x * kernelDeriv.y);
			//			v_grad2.y += SolidParticle::_mass * (v_ji.y * kernelDeriv.y);
			//		}
			//	}
			//	sf::Vector2f positionDiff = 
			//	_particles[i]._velocity = _particles[i]._velocity + (sf::Vector2f(v_grad1.x * x + v_grad2.x * y,
			//		v_grad1.y * x + v_grad2.y * y))
			//}

		}

		// Update Position and Velocity
		if (_moveParticles) { update_x_and_v(); }
	}
	else {
		EOS_solve();
	}
}

// _________________________________________________________________________________
void Simulation::run_tests() {
	if (_testNeighbors) {
		_testedParticlesId.clear();
		_testedParticlesId = TestManager::test_correct_neighbor_amount(&_particles);
	}
	if (_testKernel) {
		_testedParticlesId.clear();
		_testedParticlesId = TestManager::test_kernel(&_particles, Parameters::H);
	}
}

// _________________________________________________________________________________
void Simulation::run() {
	_totalNumSolverIterations = 0;
	_particles.reserve(_particles.size());

	_avgDensityFile.open("./avgDensityFile.dat", std::fstream::out | std::fstream::trunc);
	_timeStepFile.open("./timeStepFile.dat", std::fstream::out | std::fstream::trunc);
	_iterationsFile.open("./iterationsFile.dat", std::fstream::out | std::fstream::trunc);
	_estimatedDensityFile.open("./estimatedDensityFile.dat", std::fstream::out | std::fstream::trunc);
	if (Parameters::WRITE_SCREEN_IMAGES) {
		_renderer._screenTexture.create(Parameters::WINDOW_WIDTH, Parameters::WINDOW_HEIGHT);
		_renderer._frameCounter = 0;
	}

	if (Parameters::INTERACTIVE) {
		_videoMode = sf::VideoMode(Parameters::WINDOW_WIDTH, Parameters::WINDOW_HEIGHT);
		_window.create(_videoMode, "SPH Fluid Solver");
		sf::Time newTime;
		sf::Time renderTime;
		sf::Time lastFrameTimeStamp = _clock.getElapsedTime();
		float elapsedTime;
		float maxMaxVelocity = 0;
		float avgSolverIterations = 0;
		int frameIterations = 0;
		_maxTimeStep = 0;
		_lastTimeStep = Parameters::TIME_STEP;
		_numUpdatesPerSec = 0;
		_numIterations = 0;
		_lastSpawnTime = 0.1;
		_renderer.init_solids(&_particles);
		_watchedParticleId = (int)(_particles.size() / 2);
		while (true) {
			// Update Clock
			newTime = _clock.getElapsedTime();
			elapsedTime = newTime.asSeconds() - _lastUpdate.asSeconds();
			_numUpdatesPerSec = 1 / elapsedTime;

			_lastUpdate = _clock.getElapsedTime();

			_pauseSimulation = false;
			check_input();
			if (_endSimulation) { break; }
			if (_pauseSimulation) { continue; }

			if (_spawnParticles) {
				spawn_particles();
			}

			update_hashTable();
			_numFluidParticles = 0;
			_numParticles = _particles.size();

			_beforeNeighborhood = _clock.getElapsedTime();
			// Find Neighbors and count fluid particles
			for (int i = 0; i < _numParticles; i++) {
				if (_particles[i]._type != fluid) { continue; }
				_numFluidParticles++;
				_particles[i]._neighbors = _hashManager.return_neighbors(&_particles[i], _neighborRadius);
			}
			_neighborhoodTime += _clock.getElapsedTime() - _beforeNeighborhood;

			// Update Physics
			if (_numFluidParticles > 0) {
				_beforePhysics = _clock.getElapsedTime();
				update_physics();
				_physicsTime += _clock.getElapsedTime() - _beforePhysics;
			}

			// Update moving objects
			if (_moveSolids) { update_moving_objects(); }

			// Delete stray particles
			if (_deleteParticles) { delete_particles(); }

			if (Parameters::DOCUMENT_ITERATIONS) {
				_iterationsFile << _simulatedTime << " " << _numSolverIterations << "\n";
			}
			if (Parameters::DOCUMENT_TIME) {
				_timeStepFile << _simulatedTime << " " << _timeStepSize << "\n";
			}
			if (Parameters::DOCUMENT_AVG_DENSITY) {
				_avgDensityFile << _simulatedTime << " " << _averageDensity << "\n";
			}
			if (Parameters::DOCUMENT_ESTIMATED_DENSITY) {
				_estimatedDensityFile << _simulatedTime << " " << FluidParticle::_restDensity + _estimatedDensityError << "\n";
			}

			_simulatedTime += _timeStepSize;
			if (_simulatedTime >= Parameters::SIMULATION_LENGTH) { _endSimulation = true; }
			_cflNumber = _maxVelocity * _timeStepSize / Parameters::H;
			if (_cflNumber > _maxCflNumber) { _maxCflNumber = _cflNumber; }
			if (_cflNumber > 1) { std::cout << "CFL CONDITION VIOLATED" << std::endl; }
			if (_maxVelocity > maxMaxVelocity) { maxMaxVelocity = _maxVelocity; }
			frameIterations++;
			avgSolverIterations = _totalNumSolverIterations / frameIterations;

			if (_simulatedTime >= _nextFrame) {
				run_tests();
				_currentTime = _clock.getElapsedTime();
				if ((_currentTime - lastFrameTimeStamp).asSeconds() < _frameDistance * Parameters::SLOW_DOWN) {
					// Wait with next frame
					sf::sleep(sf::seconds(_frameDistance) * Parameters::SLOW_DOWN - (_currentTime - lastFrameTimeStamp) + sf::seconds(Parameters::TIME_OFFSET));
				}
				_renderer.update_information(_currentTime.asSeconds(),
					_simulatedTime, _particles.size(),
					_numFluidParticles, _numUpdatesPerSec, _averageDensity, _maxCflNumber,
					avgSolverIterations, _watchedParticleDensity);

				_renderer.draw(&_window, &_particles, _watchedParticleId, _markedParticlesId, _testedParticlesId);

				lastFrameTimeStamp = _currentTime;
				_nextFrame += _frameDistance;
				frameIterations = 0;
				_totalNumSolverIterations = 0;
				_maxCflNumber = 0;
			}


			// Adjust timestep if necessary
			if (Parameters::ADAPTIVE_TIME_STEP && _simulatedTime >= Parameters::INITIALIZATION_PHASE) {
				_lastTimeStep = _timeStepSize;
				_timeStepSize = Parameters::CFL_NUMBER * Parameters::H / _maxVelocity;
				if (_timeStepSize > _lastTimeStep) {
					// Cut because of last time step
					_timeStepSize -= (_timeStepSize - _lastTimeStep) * Parameters::SMOOTHING;
				}
				if (_timeStepSize > Parameters::MAX_TIME_STEP) {
					// Cut because of max time step
					_timeStepSize = Parameters::MAX_TIME_STEP;
				}
				if (_simulatedTime + _timeStepSize > _nextFrame) {
					// Cut because of next frame
					_timeStepSize = _nextFrame - _simulatedTime + Parameters::TIME_OFFSET;
				}
				else if (_simulatedTime + _timeStepSize + _timeStepSize > _nextFrame) {
					// Cut because of next 2 frames
					_timeStepSize = (_nextFrame - _simulatedTime) / 2 + Parameters::TIME_OFFSET;
				}
			}
			_numIterations++;
		}
		std::cout << "max velocity:" << maxMaxVelocity << std::endl;
		std::cout << "physics time: " << _physicsTime.asMilliseconds() << "	neighborhood time: " << _neighborhoodTime.asMilliseconds()
			<< "	total time: " << _physicsTime.asMilliseconds() + _neighborhoodTime.asMilliseconds() << std::endl;
		std::cout << "avg physics time: " << _physicsTime.asMicroseconds() / _numIterations << "	avg neighborhood time: " << _neighborhoodTime.asMicroseconds()
			/ _numIterations << "	avg total time: " << (_physicsTime.asMicroseconds() + _neighborhoodTime.asMicroseconds()) / _numIterations << std::endl;
	}

	else {
		if (!Parameters::JUST_RENDER) {
			_numIterations = 0;
			_lastSpawnTime = 0;
			_renderFile.open("./renderFile.dat", std::fstream::out | std::fstream::trunc);
			_renderFile << Parameters::TIME_STEP << std::endl;
			std::string timeString;
			std::tuple<int, int, int> rgb;
			for (int i = 0; i < _particles.size(); i++) {
				if (_particles[i]._type == solid) {
					_renderFile << "SP " << _particles[i]._position.x << " "
						<< _particles[i]._position.y << std::endl;
				}
			}

			_renderFile << "EOS" << std::endl;

			while (true) {
				_numParticles = _particles.size();

				if (_spawnParticles) {
					spawn_particles();
				}

				update_hashTable();
				_numFluidParticles = 0;
				_numParticles = _particles.size();

				// Find Neighbors and count fluid particles
				for (int i = 0; i < _numParticles; i++) {
					if (_particles[i]._type == solid) { continue; }
					_numFluidParticles++;
					_particles[i]._neighbors = _hashManager.return_neighbors(&_particles[i], _neighborRadius);
				}

				if (_numFluidParticles > 0) { update_physics(); }
				if (Parameters::DOCUMENT_ITERATIONS) {
					_iterationsFile << _simulatedTime << " " << _numSolverIterations << "\n";
				}
				if (Parameters::DOCUMENT_TIME) {
					_timeStepFile << _simulatedTime << " " << _timeStepSize << "\n";
				}
				if (Parameters::DOCUMENT_AVG_DENSITY) {
					_avgDensityFile << _simulatedTime << " " << _averageDensity << "\n";
				}
				if (Parameters::DOCUMENT_ESTIMATED_DENSITY) {
					_estimatedDensityFile << _simulatedTime << " " << FluidParticle::_restDensity + _estimatedDensityError << "\n";
				}

				_simulatedTime += _timeStepSize;
				_cflNumber = _maxVelocity * _timeStepSize / Parameters::H;

				if (_simulatedTime >= _nextFrame) {
					timeString = std::to_string(_simulatedTime);
					timeString.resize(5, ' ');
					std::cout << "Calculated time: " << timeString << " / "
						<< Parameters::SIMULATION_LENGTH << "\r";
					for (int i = 0; i < _numParticles; i++) {
						if (_particles[i]._type == fluid && Parameters::COLOR_CODE_SPEED) {
							_renderFile << "FP " << _particles[i]._position.x << " " << _particles[i]._position.y << " " <<
								0 << " " << _particles[i]._colorFactor << " " << 255 << std::endl;
						}
						else if (_particles[i]._type == fluid && Parameters::COLOR_CODE_PRESSURE) {
							rgb = Functions::color_code_pressure(_particles[i]._pressure);
							_renderFile << "FP " << _particles[i]._position.x << " " << _particles[i]._position.y << " " <<
								std::get<0>(rgb) << " " << std::get<1>(rgb) << " " << std::get<2>(rgb) << std::endl;
						}
						else if (_particles[i]._type == fluid) {
							_renderFile << "FP " << _particles[i]._position.x << " " << _particles[i]._position.y << " " <<
								0 << " " << 0 << " " << 255 << std::endl;
						}
					}
					_renderFile << "D " << _averageDensity << "\n";
					_renderFile << "EOU" << std::endl;
					_nextFrame += _frameDistance;
				}
	
				// Adjust timestep if necessary
				if (Parameters::ADAPTIVE_TIME_STEP && _simulatedTime > Parameters::INITIALIZATION_PHASE) {
					_timeStepSize = Parameters::CFL_NUMBER * Parameters::H / _maxVelocity;
					if (_timeStepSize > Parameters::MAX_TIME_STEP) {
						_timeStepSize = Parameters::MAX_TIME_STEP;
					}
					if (_simulatedTime + _timeStepSize > _nextFrame) {
						_timeStepSize = _nextFrame - _simulatedTime + Parameters::TIME_OFFSET;
					}
					else if (_simulatedTime + _timeStepSize + _timeStepSize > _nextFrame) {
						_timeStepSize = (_nextFrame - _simulatedTime + Parameters::TIME_OFFSET) / 2;
					}
				}
				_numIterations++;
				if (_simulatedTime > Parameters::SIMULATION_LENGTH) { break; }
			}
			_renderFile << "Z" << std::endl;
			_renderFile.close();
			_avgDensityFile.close();
			_timeStepFile.close();
			_iterationsFile.close();
			_estimatedDensityFile.close();
			std::cout << _avgNeighborhoodTime / Parameters::SIMULATION_LENGTH << std::endl;
		}

		_renderFile.open("./renderFile.dat", std::fstream::in);
		_videoMode = sf::VideoMode(Parameters::WINDOW_WIDTH, Parameters::WINDOW_HEIGHT);


		sf::Time elapsedTime;
		sf::Time timeDifference;
		int numUpdatesPerSec;
		int numShapes = 0;
		int numFluids = 0;
		int numSolids = 0;
		std::string type;
		float xValue;
		float yValue;
		int r;
		int g;
		int b;
		float density;
		float timeStepSize;
		int c = 0;
		_endSimulation = false;
		_lastUpdate = sf::seconds(0);
		_renderFile >> timeStepSize;
		_clock.restart();

		std::vector<sf::CircleShape> solidShapes = std::vector<sf::CircleShape>();
		std::vector<std::vector<sf::CircleShape>> fluidShapes = std::vector<std::vector<sf::CircleShape>>();
		std::vector<sf::CircleShape> currentFluids = std::vector<sf::CircleShape>();
		std::vector<float> densities = std::vector<float>();
		sf::CircleShape newShape = sf::CircleShape();

		// Solid reading Loop
		std::cout << "reading solids\r";
		while (true) {
			_renderFile >> type;
			if (type[0] == 'E') { break; }
			if (type[0] == 'Z') { _endSimulation = true; break; }
			_renderFile >> xValue >> yValue;
			newShape.setFillColor(sf::Color::Color::White);
			newShape.setRadius(SolidParticle::_size * _zoomFactor);
			newShape.setPosition(xValue * _zoomFactor, yValue * _zoomFactor);
			solidShapes.push_back(newShape);
			numSolids++;
		}

		// Fluid Reading Loop
		std::cout << std::endl;
		while(true) {
			if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
				// Close application
				_endSimulation = true;
				break;
			}
			_renderFile >> type;
			if (type[0] == 'Z') { break; }
			if (type[0] == 'E') {
				fluidShapes.push_back(currentFluids);
				int size = currentFluids.size();
				currentFluids.clear();
				currentFluids.reserve(size);
				c++;
				std::cout << "Reading Frame: " << c << " of " << Parameters::SIMULATION_LENGTH / _frameDistance << "\r";
				continue;
			}
			if (type[0] == 'D') {
				_renderFile >> density;
				densities.push_back(density);
				continue;
			}
			_renderFile >> xValue >> yValue >> r >> g >> b;

			if (type[0] == 'F') {
				newShape.setRadius(FluidParticle::_size * _zoomFactor / 2);
				newShape.setFillColor(sf::Color::Color(r, g, b));
				newShape.setPosition(xValue* _zoomFactor, yValue* _zoomFactor);
				currentFluids.push_back(newShape);
			}
			_renderFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		std::cout << "=======================================================" << std::endl;
		std::cout << "COMPUTATIONS READY! Press Space to watch the Simulation" << std::endl;
		std::cout << "=======================================================" << std::endl;

		while (!sf::Keyboard::isKeyPressed(sf::Keyboard::Space)) {
			continue;
		}
		_window.create(_videoMode, "SPH Fluid Solver");
		while (!_endSimulation) {

			_clock.restart();
			_lastUpdate = sf::seconds(0);
			// Render Loop
			for (int i = 0; i < c; i++) {
				if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
					// Close application
					_endSimulation = true;
					break;
				}
				elapsedTime = _clock.getElapsedTime();
				timeDifference = elapsedTime - _lastUpdate;
				if (timeDifference.asSeconds() < _frameDistance * Parameters::SLOW_DOWN) {
					// Wait with next frame
					sf::sleep(sf::seconds(_frameDistance) * Parameters::SLOW_DOWN - timeDifference + sf::seconds(Parameters::TIME_OFFSET));
				}
				numFluids = fluidShapes[i].size();
				numShapes = numFluids + numSolids;
				_averageDensity = densities[i];
				numUpdatesPerSec = 1 / (elapsedTime.asSeconds() - _lastUpdate.asSeconds());
				_lastUpdate = _clock.getElapsedTime();
				_window.clear();

				for (int ii = 0; ii < numFluids; ii++) {
					_window.draw(fluidShapes[i][ii]);
				}
				for (int ii = 0; ii < numSolids; ii++) {
					_window.draw(solidShapes[ii]);
				}

				_renderer.update_information(elapsedTime.asSeconds(), i * _frameDistance,
					numShapes, numFluids, numUpdatesPerSec, _averageDensity,
					Parameters::CFL_NUMBER);
				_window.draw(_renderer._infoPanel);
				_window.draw(_renderer._description);
				_window.draw(_renderer._information);
				_window.display();
			}

		}
		_window.close();
	}

	_renderFile.close();
	_timeStepFile.close();
	_iterationsFile.close();
	_avgDensityFile.close();
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

	for (int timeSteps = 0; timeSteps < Parameters::SIMULATION_LENGTH; timeSteps++) {
		
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
				_avgDensityFile << c * Parameters::TIME_STEP << " " << _averageDensity << "\n";
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
	_timeStepFile.close();
	_iterationsFile.close();
	_avgDensityFile.close();

}