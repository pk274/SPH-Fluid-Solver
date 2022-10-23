// Paul Kull, 2022

#include <cstdlib>
#include <ctime>

#include "./Simulation.h"

// _________________________________________________________________________________
Simulation::Simulation(SimulationPreset preset = StuffedBox, int framelimit) {

	_particles = std::vector<Particle>();
	_particles.clear();

	_neighborRadius = 4;

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

// __________________________________________________________________________________
float calculate_distance(sf::Vector2f pos1, sf::Vector2f pos2) {
	return std::abs(std::sqrt(std::pow(pos1.x - pos2.x, 2) + std::pow(pos1.y - pos2.y, 2)));
}


// _________________________________________________________________________________
std::vector<Particle*> Simulation::find_neighbors(Particle* particle) {
	std::vector<Particle*> neighbors = std::vector<Particle*>();
	std::vector<Particle*> possibleNeighbors = _hashManager.return_possible_neighbors(particle);
	
	int numPosNeighs = possibleNeighbors.size();
	for (int i = 0; i < numPosNeighs; i++) {
		if (possibleNeighbors[i]->_id == particle->_id) { goto END_OF_OUTER_LOOP; }
		for (int j = 0; j < neighbors.size(); j++) {
			if (possibleNeighbors[i]->_id == neighbors[j]->_id) { goto END_OF_OUTER_LOOP; }
		}
		if (calculate_distance(possibleNeighbors[i]->_position, particle->_position) <= _neighborRadius) {
			neighbors.push_back(possibleNeighbors[i]);
		}
	END_OF_OUTER_LOOP:;
	}
	return neighbors;

}

// _________________________________________________________________________________
void Simulation::update_physics() {
	_markedParticlesId.clear();
	int numParticles = _particles.size();
	for (int i = 0; i < numParticles; i++) {
		std::vector<Particle*> neighbors = find_neighbors(&_particles[i]);
		if (_particles[i]._id == _watchedParticleId) {
			for (int j = 0; j < neighbors.size(); j++) {
				_markedParticlesId.push_back(neighbors[j]->_id);
			}
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