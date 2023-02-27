// Paul Kull, 2023

#include "./InitManager.h"


// Helpful Structures to place

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
std::vector<Particle> placePool(sf::Vector2f pos, int size, int startId = 0) {
	std::vector<Particle> box = std::vector<Particle>();
	for (int i = 0; i < size / 2; i++) {
		box.push_back(SolidParticle(startId + box.size(), pos));
		pos.y += SolidParticle::_size;
	}
	for (int i = 0; i < size; i++) {
		box.push_back(SolidParticle(startId + box.size(), pos));
		pos.x += SolidParticle::_size;
	}
	for (int i = 0; i < size / 2; i++) {
		box.push_back(SolidParticle(startId + box.size(), pos));
		pos.y -= SolidParticle::_size;
	}
	return box;
}

// _________________________________________________________________________________
std::vector<Particle> placeParticleLine(sf::Vector2f pos1, sf::Vector2f pos2, ParticleType type,
	int startId = 0, bool placeFirst = 1, bool placeLast = 1) {

	std::vector<Particle> line = std::vector<Particle>();
	sf::Vector2f distance = Functions::calculate_distance(pos1, pos2);
	float distanceNorm = Functions::calculate_distance_norm(distance);
	float particlesNecessary = distanceNorm;
	sf::Vector2f step;
	if (type == solid) {
		particlesNecessary = distanceNorm / SolidParticle::_size;
		step = -distance * (SolidParticle::_size / distanceNorm);
	}
	else if (type == fluid) {
		particlesNecessary = distanceNorm / FluidParticle::_size;
		step = -distance * (FluidParticle::_size / distanceNorm);
	}
	sf::Vector2f pos = pos1;

	for (int i = 0; i < particlesNecessary; i++) {
		if (i == 0 && !placeFirst) { pos += step; continue; }
		if (type == solid) {
			line.push_back(SolidParticle(startId + line.size(), pos));
		}
		else if (type == fluid) {
			line.push_back(FluidParticle(startId + line.size(), pos));
		}
		pos += step;
	}
	if (!placeLast) { line.pop_back(); }
	return line;
}

// _________________________________________________________________________________
std::vector<Particle> placeTriangle(sf::Vector2f pos1, sf::Vector2f pos2, sf::Vector2f pos3, ParticleType type,
	int startId = 0, bool placeFirst = 1, bool placeLast = 1) {
	std::vector<Particle> triangle = std::vector<Particle>();
	std::vector<Particle> line;
	if (type == solid) {
		line = placeParticleLine(pos1, pos2, solid, startId + triangle.size(), 1, 1);
		for (int i = 0; i < line.size(); i++) {
			triangle.push_back(line[i]);
		}
		line = placeParticleLine(pos2, pos3, solid, startId + triangle.size(), 1, 1);
		for (int i = 0; i < line.size(); i++) {
			triangle.push_back(line[i]);
		}
		line = placeParticleLine(pos3, pos1, solid, startId + triangle.size(), 1, 1);
		for (int i = 0; i < line.size(); i++) {
			triangle.push_back(line[i]);
		}
	}
	return triangle;
}


// Initialization procedures:

// _________________________________________________________________________________
InitManager::InitManager(Simulation* sim) {
	_sim = sim;
}

// _________________________________________________________________________________
void InitManager::init_simulation(SimulationPreset preset) {
	
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
		init_breaking_dam_simulation(80, 5);
		break;
	case BigBreakingDam:
		init_breaking_dam_simulation(150, 3);
		break;
	case GiantBox:
		init_stuffed_box_simulation(100, 3);
		break;
	case FourLayers:
		init_layer_simulation(39, 10, 4);
		break;
	case ManyLayers:
		init_layer_simulation(60, 7, 10);
		break;
	case WideLayers:
		init_wide_layer_simulation(150, 2, 10);
		break;
	case Cup:
		init_cup_simulation(120, 4);	//3
		break;
	case Complex:
		init_complex_simulation(200, 2);
		break;
	case Osmosis:
		init_osmosis_simulation(100, 4);
		break;
	}
}



// _________________________________________________________________________________
void InitManager::init_empty_simulation() {

	_sim->_videoMode = sf::VideoMode(Parameters::WINDOW_WIDTH, Parameters::WINDOW_HEIGHT);
	_sim->_zoomFactor = 50;
	_sim->_window.create(_sim->_videoMode, "SPH Fluid Solver");
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);

	sf::Vector2f pos = sf::Vector2f(5, 5);
	_sim->_particles.push_back(Particle(fluid, 0, pos));
	_sim->_moveParticles = false;
}


// _________________________________________________________________________________
void  InitManager::init_stuffed_box_simulation(int size, int zoom) {

	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);
	// _hashManager = HashManager(_neighborRadius, size / 2, size / 2);

	// Add Particles
	// Add Particles for arena
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeBox(pos, size - 2, _sim->_particles.size());
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * 2);
	box = placeBox(pos, size - 4, _sim->_particles.size());
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}

	pos = sf::Vector2f(SolidParticle::_size * 3, SolidParticle::_size * 3);

	for (int i = 0; i < size - 5; i++) {
		for (int ii = 0; ii < size - 5; ii++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.y += FluidParticle::_size;
		}
		pos.y = SolidParticle::_size * 3;
		pos.x += FluidParticle::_size;
	}
	_sim->_moveParticles = false;
	_sim->_testNeighbors = true;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = false;
}

// ______________________________________________________________________________________________________
void InitManager::init_single_particle_simulation(int size, int zoom) {

	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);
	// _hashManager = HashManager(_sim->_neighborRadius, size / 2, size / 2);

	// Add Particles
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeBox(pos, size - 2, _sim->_particles.size());
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}


	pos.x = 16;
	pos.y = 35;
	_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));

	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = true;
	_sim->_deleteParticles = true;
}

// _________________________________________________________________________________
void InitManager::init_rotated_box_simulation(int size, int zoom, int rotation) {

	init_stuffed_box_simulation(size, zoom);
	sf::Vector2f offset = sf::Vector2f(0, size * 1.5);
	for (int i = 0; i < _sim->_particles.size(); i++) {
		_sim->_particles[i]._position.x += std::cos(rotation * 3.141592653589 / 180) * _sim->_particles[i]._position.x
			+ std::sin(rotation * 3.141592653589 / 180) * _sim->_particles[i]._position.y;
		_sim->_particles[i]._position.y += -std::sin(rotation * 3.141592653589 / 180) * _sim->_particles[i]._position.x
			+ std::cos(rotation * 3.141592653589 / 180) * _sim->_particles[i]._position.y;
		_sim->_particles[i]._position += offset;
	}
	_sim->_testNeighbors = true;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = false;
}

// _________________________________________________________________________________
void InitManager::init_random_particles_simulation(int size, int zoom, int numParticles) {

	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);
	// _hashManager = HashManager(_sim->_neighborRadius, size / 2, size / 2);

	// Add Particles for arena
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeBox(pos, size - 2, _sim->_particles.size());
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	for (int i = 0; i < numParticles; i++) {
		float x = (std::rand() % (size - 10)) + SolidParticle::_size * 5;
		float y = (std::rand() % (size - 10)) + SolidParticle::_size * 5;
		_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), sf::Vector2f(x, y)));
	}
	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = false;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
}


// _________________________________________________________________________________
void  InitManager::init_breaking_dam_simulation(int size, int zoom) {

	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);
	// _hashManager = HashManager(_neighborRadius, size / 2, size / 2);

	// Add Particles for arena
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeBox(pos, size - 2, _sim->_particles.size());
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * (size - 2));
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		pos.x = SolidParticle::_size * 2;
		pos.y -= FluidParticle::_size;
	}
	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;

}

// _________________________________________________________________________________
void InitManager::init_layer_simulation(int size, int zoom, int layers) {
	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 4, size / 2);
	// _hashManager = HashManager(_neighborRadius, size / 4, size / 2);

	// Add Particles for arena
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeTallBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeTallBox(pos, size - 2, 1, _sim->_particles.size());
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * (size - 2));
	for (int i = 0; i < layers; i++) {
		for (int j = 0; j < size / 2 - 3; j++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		pos.x = SolidParticle::_size * 2;
		pos.y -= FluidParticle::_size;
	}
	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = false;
}




// _________________________________________________________________________________
void InitManager::init_wide_layer_simulation(int size, int zoom, int layers) {
	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);
	// _hashManager = HashManager(_neighborRadius, size / 2, size / 2);

	// Add Particles for arena
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> box = placeBox(pos, size);
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeBox(pos, size - 2, 1);
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * (size - 2));
	for (int i = 0; i < layers; i++) {
		for (int j = 0; j < size - 3; j++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		pos.x = SolidParticle::_size * 2;
		pos.y -= FluidParticle::_size;
	}
	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = false;
}


// _________________________________________________________________________________
void InitManager::init_cup_simulation(int size, int zoom) {
	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);
	//_hashManager = HashManager(_neighborRadius, size / 2, size / 2);

	int downshift = 130;
	int sideshift = 12;

	// Add Particles for arena
	sf::Vector2f pos1 = sf::Vector2f(sideshift, downshift);
	sf::Vector2f pos2 = sf::Vector2f(sideshift, downshift);
	std::vector<Particle> box = placePool(pos1, size);
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos1 = sf::Vector2f(SolidParticle::_size + sideshift, downshift);
	box = placePool(pos1, size - 2, _sim->_particles.size());
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos1 = sf::Vector2f(50, 10);
	std::vector<Particle> line1 = placeParticleLine(pos1, pos1 + sf::Vector2f(90, 50), solid, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}
	pos2 = pos1 + sf::Vector2f(-25, 45);
	line1 = placeParticleLine(pos2, pos2 + sf::Vector2f(90, 50), solid, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}

	std::vector<Particle> line2;
	line1 = placeParticleLine(pos2, pos1, solid, _sim->_particles.size(), 0, 1);
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}
	for (int i = 0; i < line1.size(); i++) {
		line2 = placeParticleLine(line1[i]._position, line1[i]._position + sf::Vector2f(90, 50), fluid, _sim->_particles.size(), 0, 0);
		for (int j = 0; j < line2.size(); j++) {
			_sim->_particles.push_back(line2[j]);
		}
	}

	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = false;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
}

// _________________________________________________________________________________
void InitManager::init_complex_simulation(int size, int zoom) {
	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);
	//_hashManager = HashManager(_neighborRadius, size / 2, size / 2);

	// Add Particles for arena
	sf::Vector2f pos1 = sf::Vector2f(0, 0);
	sf::Vector2f pos2 = sf::Vector2f(0, 0);
	sf::Vector2f pos3 = sf::Vector2f(0, 0);
	std::vector<Particle> shape;
	std::vector<Particle> box = placeBox(pos1, size);
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}
	pos1 = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeBox(pos1, size - 2, _sim->_particles.size());
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}

	pos1 = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * 2);
	for (int i = 0; i < size / 3; i++) {
		for (int j = 0; j < size / 3; j++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos1));
			pos1.x += FluidParticle::_size;
		}
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos1));
		pos3 = pos1;
		pos1.x = SolidParticle::_size * 2;
		pos1.y += FluidParticle::_size;
	}
	pos2 = pos1 + sf::Vector2f(SolidParticle::_size * size / 1.5, SolidParticle::_size * size / 4);
	shape = placeParticleLine(pos1, pos2, solid, _sim->_particles.size());
	for (int i = 0; i < shape.size(); i++) {
		_sim->_particles.push_back(shape[i]);
	}
	pos2 += sf::Vector2f(0, -SolidParticle::_size * 8);
	shape = placeParticleLine(pos3, pos2, solid, _sim->_particles.size());
	for (int i = 0; i < shape.size(); i++) {
		_sim->_particles.push_back(shape[i]);
	}

	pos1 = sf::Vector2f(SolidParticle::_size * (size - SolidParticle::_size * 5), SolidParticle::_size * size / 2);
	pos2 = sf::Vector2f(SolidParticle::_size * (size - SolidParticle::_size * 5), SolidParticle::_size * size - 10);
	pos3 = sf::Vector2f(SolidParticle::_size * (size - SolidParticle::_size * 10), SolidParticle::_size * size - 10);
	shape = placeTriangle(pos1, pos2, pos3, solid, _sim->_particles.size());
	for (int i = 0; i < shape.size(); i++) {
		_sim->_particles.push_back(shape[i]);
	}

	pos1 = sf::Vector2f(SolidParticle::_size * (size / 6), SolidParticle::_size * size * 3 / 4);
	pos2 = sf::Vector2f(SolidParticle::_size * (size / 6), SolidParticle::_size * size - 10);
	pos3 = sf::Vector2f(SolidParticle::_size * (size / 3), SolidParticle::_size * size - 10);
	shape = placeTriangle(pos1, pos2, pos3, solid, _sim->_particles.size());
	for (int i = 0; i < shape.size(); i++) {
		_sim->_particles.push_back(shape[i]);
	}

	pos1 = sf::Vector2f(SolidParticle::_size * (size / 2 + 10), SolidParticle::_size * size * 2 / 3);
	pos2 = sf::Vector2f(SolidParticle::_size * (size / 2 - 10), SolidParticle::_size * size * 2 / 3);
	pos3 = sf::Vector2f(SolidParticle::_size * (size / 2), SolidParticle::_size * size * 2 / 2.5);
	shape = placeTriangle(pos1, pos2, pos3, solid, _sim->_particles.size());
	for (int i = 0; i < shape.size(); i++) {
		_sim->_particles.push_back(shape[i]);
	}


	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
}

// _________________________________________________________________________________
void InitManager::init_osmosis_simulation(int size, int zoom) {
	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);
	// _hashManager = HashManager(_neighborRadius, size / 2, size / 2);

	// Add Particles for arena
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> shape = placeBox(pos, size);
	for (int i = 0; i < shape.size(); i++) {
		_sim->_particles.push_back(shape[i]);
	}
	pos = sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	shape = placeBox(pos, size - 2, _sim->_particles.size());
	for (int i = 0; i < shape.size(); i++) {
		_sim->_particles.push_back(shape[i]);
	}

	pos = sf::Vector2f(SolidParticle::_size * size / 2 + 1, SolidParticle::_size * (size - 5));
	shape = placeParticleLine(pos, pos - sf::Vector2f(0, size * 2 - 6), solid, _sim->_particles.size());
	for (int i = 0; i < shape.size(); i++) {
		_sim->_particles.push_back(shape[i]);
	}

	pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * (size - 2));
	for (int i = 0; i < size / 1.5; i++) {
		for (int j = 0; j < size / 2 - 1; j++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		pos.x = SolidParticle::_size * 2;
		pos.y -= FluidParticle::_size;
	}

	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
}