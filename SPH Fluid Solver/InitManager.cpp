// Paul Kull, 2023

#include <iostream>
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
	if (type == solid || type == moving) {
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
		else if (type == moving) {
			line.push_back(SolidParticle(startId + line.size(), pos));
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

// _________________________________________________________________________________
std::vector<Particle> place_fountain(sf::Vector2f lowerLeft, Simulation* sim, int startId = 0) {
	int height = 400;
	int width = 350;
	sf::Vector2f pos1 = lowerLeft;
	sf::Vector2f pos2 = lowerLeft;
	sf::Vector2f pos3 = lowerLeft;
	sf::Vector2f pos4 = lowerLeft;
	std::vector<Particle> object1 = std::vector<Particle>();
	std::vector<Particle> object2 = std::vector<Particle>();
	std::vector<Particle> fountain = std::vector<Particle>();

	// Place first pond
	pos1.y -= height * SolidParticle::_size / 10;
	for (int i = 0; i < height / 10; i++) {
		fountain.push_back(SolidParticle(startId + fountain.size(), pos1));
		pos1.y += SolidParticle::_size;
	}
	for (int i = 0; i < width; i++) {
		fountain.push_back(SolidParticle(startId + fountain.size(), pos1));
		pos1.x += SolidParticle::_size;
	}
	for (int i = 0; i < height / 10; i++) {
		fountain.push_back(SolidParticle(startId + fountain.size(), pos1));
		pos1.y -= SolidParticle::_size;
	}
	pos1 = lowerLeft;
	pos1.y -= height * SolidParticle::_size / 10;
	pos1.x += SolidParticle::_size;
	for (int i = 0; i < height / 10 - 1; i++) {
		fountain.push_back(SolidParticle(startId + fountain.size(), pos1));
		pos1.y += SolidParticle::_size;
	}
	for (int i = 0; i < width - 2; i++) {
		fountain.push_back(SolidParticle(startId + fountain.size(), pos1));
		pos1.x += SolidParticle::_size;
	}
	for (int i = 0; i < height / 10 - 1; i++) {
		fountain.push_back(FluidParticle(startId + fountain.size(), pos1));
		pos1.y -= SolidParticle::_size;
	}

	// Second pond
	pos1.x = lowerLeft.x + width * SolidParticle::_size / 5;
	pos1.y = lowerLeft.y - height * SolidParticle::_size * 7/20;
	pos2.x = lowerLeft.x + width * SolidParticle::_size / 4;
	pos2.y = lowerLeft.y - height * SolidParticle::_size * 7 / 20 + SolidParticle::_size * height / 12;
	object1 = placeParticleLine(pos1, pos2, solid, fountain.size(), true, true);
	for (int i = 0; i < object1.size(); i++) {
		fountain.push_back(object1[i]);
	}
	for (int i = 0; i < width / 2; i++) {
		fountain.push_back(SolidParticle(startId + fountain.size(), pos2));
		pos2.x += SolidParticle::_size;
	}
	pos1.x += 6 * width / 5;
	object1 = placeParticleLine(pos1, pos2, solid, fountain.size(), true, true);
	for (int i = 0; i < object1.size(); i++) {
		fountain.push_back(object1[i]);
	}

	// Third pond
	pos1.x = lowerLeft.x + width * SolidParticle::_size / 3;
	pos1.y = lowerLeft.y - height * SolidParticle::_size * 6/10;
	pos2.x = lowerLeft.x + width * SolidParticle::_size * 23 / 60;
	pos2.y = lowerLeft.y - height * SolidParticle::_size * 6 / 10 + height * SolidParticle::_size / 15;
	object1 = placeParticleLine(pos1, pos2, solid, fountain.size(), true, true);
	for (int i = 0; i < object1.size(); i++) {
		fountain.push_back(object1[i]);
	}
	for (int i = 0; i < width * 7 / 30; i++) {
		fountain.push_back(SolidParticle(startId + fountain.size(), pos2));
		pos2.x += SolidParticle::_size;
	}
	pos1.x += 2 * width / 3;
	object1 = placeParticleLine(pos1, pos2, solid, fountain.size(), true, true);
	for (int i = 0; i < object1.size(); i++) {
		fountain.push_back(object1[i]);
	}

	// Place central pillar
	pos1.x = lowerLeft.x + width * SolidParticle::_size / 2 - width * SolidParticle::_size / 100;
	pos1.y = lowerLeft.y;
	for (int i = 0; i < height * 11/16; i++) {
		fountain.push_back(SolidParticle(startId + fountain.size(), pos1));
		pos1.y -= SolidParticle::_size;
	}
	pos1.x = lowerLeft.x + width * SolidParticle::_size / 2 + width * SolidParticle::_size / 100;
	pos1.y = lowerLeft.y;
	for (int i = 0; i < height * 11/16; i++) {
		fountain.push_back(SolidParticle(startId + fountain.size(), pos1));
		pos1.y -= SolidParticle::_size;
	}

	int openPillar = -SolidParticle::_size * 30;
	pos1.x = lowerLeft.x + width * SolidParticle::_size / 2 - width * SolidParticle::_size / 100;
	pos1.y = lowerLeft.y - (height * SolidParticle::_size * 11 / 16) + 2 * SolidParticle::_size - openPillar;
	pos1.x += SolidParticle::_size;
	for (int i = 0; i < (width / 50) - 1; i++) {
		fountain.push_back(SolidParticle(startId + fountain.size(), pos1));
		pos1.x += SolidParticle::_size;
	}

	// Add central Spawns
	pos1.x = lowerLeft.x + width * SolidParticle::_size / 2 - width * SolidParticle::_size / 100;
	pos1.y = lowerLeft.y - (height * SolidParticle::_size * 11 / 16) + SolidParticle::_size - openPillar;
	pos1.x += SolidParticle::_size;
	pos2.x = lowerLeft.x + width * SolidParticle::_size / 2;
	pos2.y = lowerLeft.y - (height * SolidParticle::_size * 11 / 16) + SolidParticle::_size - openPillar;
	for (int i = 0; i < (width / 50) - 1; i++) {
		sim->_spawnLocations.push_back(pos1);
		sim->_spawnVelocities.push_back(sf::Vector2f(0, - height * 1.2
			- (height / 5 / (1.5 + Functions::calculate_distance_norm(Functions::calculate_distance(pos1, pos2))))));
		pos1.x += SolidParticle::_size;
	}

	// Add side spawns
	// left
	pos1 = lowerLeft;
	pos1.y -= height * SolidParticle::_size / 10 + 1;
	sim->_spawnLocations.push_back(pos1);
	sim->_spawnVelocities.push_back(sf::Vector2f(width / 2.5, -height * 2.2));
	pos1.x += SolidParticle::_size;
	sim->_spawnLocations.push_back(pos1);
	sim->_spawnVelocities.push_back(sf::Vector2f(width / 2.5, -height * 2.2));
	// right
	pos1 = lowerLeft;
	pos1.x += width * SolidParticle::_size;
	pos1.y -= height * SolidParticle::_size / 10 + 1;
	sim->_spawnLocations.push_back(pos1);
	sim->_spawnVelocities.push_back(sf::Vector2f(- width / 2.5, -height * 2.2));
	pos1.x -= SolidParticle::_size;
	sim->_spawnLocations.push_back(pos1);
	sim->_spawnVelocities.push_back(sf::Vector2f(- width / 2.5, -height * 2.2));


	sim->_spawnDelay = 0.005;


	return fountain;
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
		init_breaking_dam_simulation(60, 8);
		break;
	case BigBreakingDam:
		init_breaking_dam_simulation(120, 4);
		break;
	case FixedBreakingDam:
		init_breaking_dam_simulation(200, 2, false, 10000);
		break;
	case SmallBreakingDam:
		init_breaking_dam_simulation(40, 13);
		break;
	case TallBreakingDam:
		init_tall_breaking_dam_simulation(250, 2);
		break;
	case GiantBox:
		init_stuffed_box_simulation(50, 4);
		break;
	case FourLayers:
		init_layer_simulation(39, 10, 4);
		break;
	case ManyLayers:
		init_layer_simulation(59, 8, 40, 1, 1);
		break;
	case WideLayers:
		init_wide_layer_simulation(80, 6, 40);	// 50 for error analysis, 80 for layers
		break;
	case Cup:
		init_cup_simulation(120, 4);	//3
		break;
	case Complex:
		init_complex_simulation(240, 2);
		break;
	case Complex2:
		init_complex_2_simulation(240, 2);
		break;
	case Osmosis:
		init_osmosis_simulation(100, 4);
		break;
	case SideSpawn:
		init_side_spawn_simulation(100, 4);
		break;
	case Rain:
		init_rain_simulation(70, 5);
		break;
	case Rain2:
		init_rain2_simulation(70, 5);
		break;
	case Fountain:
		init_fountain_simulation();
		break;
	case MovingWall:
		init_moving_wall_simulation(50, 8);
		break;
	case WaveGenerator:
		init_wave_generator_simulation(120, 4);
		break;
	case CubeDrop:
		init_cube_drop_simulation(120, 4);
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
	_sim->_spawnParticles = false;
}


// _________________________________________________________________________________
void  InitManager::init_stuffed_box_simulation(int size, int zoom) {

	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);
	// _hashManager = HashManager(_neighborRadius, size / 2, size / 2);

	// Add Particles
	// Add Particles for arena
	sf::Vector2f pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * 2);
	std::vector<Particle> box = placeBox(pos, size - 4, _sim->_particles.size());
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
	}

	pos = sf::Vector2f(SolidParticle::_size * 3, SolidParticle::_size * 3);

	for (int i = 0; i < size - 5; i++) {
		for (int ii = 0; ii < size - 5; ii++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		pos.x = SolidParticle::_size * 3;
		pos.y += FluidParticle::_size;
	}
	_sim->_moveParticles = false;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = true;
	_sim->_deleteParticles = false;
	_sim->_spawnParticles = false;
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
	pos.y = size - 3.0;
	_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos, sf::Vector2f(8, 0)));

	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = true;
	_sim->_deleteParticles = true;
	_sim->_spawnParticles = false;
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
	_sim->_spawnParticles = false;
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
		_sim->_particles.at(_sim->_particles.size() - 1)._velocity.x = 10.f;
		//_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), sf::Vector2f(x + 3, y)));
		//_sim->_particles.at(_sim->_particles.size() - 1)._velocity.x = -10.f;
	}
	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = false;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
	_sim->_spawnParticles = false;
}


// _________________________________________________________________________________
void  InitManager::init_breaking_dam_simulation(int size, int zoom, bool adaptive, int numFluidParticles) {

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

	float xOffsetSize = 0.02;	// 0.3
	float yOffsetSize = 0.02;	// 0.05
	float xFluidSize = 1;
	float yFluidSize = 1;
	pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * (size - 2));
	pos.x = SolidParticle::_size * 2 + xOffsetSize;
	if (adaptive) { xFluidSize = size / 2.2; yFluidSize = size / 1.8; }
	else { xFluidSize = std::sqrt(numFluidParticles); yFluidSize = std::sqrt(numFluidParticles); }
	for (int i = 0; i < yFluidSize; i++) {
		for (int j = 0; j < xFluidSize; j++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		xOffsetSize = xOffsetSize * -1;
		pos.x = SolidParticle::_size * 2 + xOffsetSize;
		pos.y -= FluidParticle::_size + yOffsetSize;
	}
	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
	_sim->_spawnParticles = false;

}

// _________________________________________________________________________________
void  InitManager::init_tall_breaking_dam_simulation(int size, int zoom) {

	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 4, size / 2);
	// _hashManager = HashManager(_neighborRadius, size / 2, size / 2);

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

	float xOffsetSize = 0.02;	// 0.3
	float yOffsetSize = 0.02;	// 0.05
	float xFluidSize = 1;
	float yFluidSize = 1;
	pos = sf::Vector2f(SolidParticle::_size * 2, SolidParticle::_size * (size - 2));
	pos.x = SolidParticle::_size * 2 + xOffsetSize;
	xFluidSize = size / 7; yFluidSize = size / 1.5;
	for (int i = 0; i < yFluidSize; i++) {
		for (int j = 0; j < xFluidSize; j++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		xOffsetSize = xOffsetSize * -1;
		pos.x = SolidParticle::_size * 2 + xOffsetSize;
		pos.y -= FluidParticle::_size + yOffsetSize;
	}
	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
	_sim->_spawnParticles = false;

}

// _________________________________________________________________________________
void InitManager::init_layer_simulation(int size, int zoom, int layers, bool xOffset, bool yOffset) {
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
	float xOffsetSize = 0.5;
	float yOffsetSize = 0.09;	//0.05
	if (!xOffset) { xOffsetSize = 0; }
	if (!yOffset) { yOffsetSize = 0; }
	pos = sf::Vector2f(SolidParticle::_size * 2 + xOffsetSize / 2, SolidParticle::_size * (size - 2));
	for (int i = 0; i < layers; i++) {
		for (int j = 0; j < size / 2 - 3; j++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		xOffsetSize = xOffsetSize * -1;
		pos.x = SolidParticle::_size * 2 + xOffsetSize;
		pos.y -= FluidParticle::_size + yOffsetSize;
	}
	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = false;
	_sim->_spawnParticles = false;
}




// _________________________________________________________________________________
void InitManager::init_wide_layer_simulation(int size, int zoom, int layers, bool xOffset, bool yOffset) {
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
	float xOffsetSize = 0.5;
	float yOffsetSize = 0.05;
	if (!xOffset) { xOffsetSize = 0; }
	if (!yOffset) { yOffsetSize = 0; }
	pos = sf::Vector2f(SolidParticle::_size * 2 + xOffsetSize / 2, SolidParticle::_size * (size - 2));
	for (int i = 0; i < layers; i++) {
		for (int j = 0; j < size - 3; j++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		xOffsetSize = xOffsetSize * -1;
		pos.x = SolidParticle::_size * 2 + xOffsetSize;
		pos.y -= FluidParticle::_size + yOffsetSize;
	}
	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = false;
	_sim->_spawnParticles = false;
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
	_sim->_spawnParticles = false;
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
	_sim->_spawnParticles = false;
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
		for (int j = 0; j < size / 2 - 2; j++) {
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
	_sim->_spawnParticles = false;
}


// _________________________________________________________________________________
void InitManager::init_side_spawn_simulation(int size, int zoom) {

	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);

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

	_sim->_spawnLocations.push_back(sf::Vector2f(SolidParticle::_size * 4, SolidParticle::_size * 5));
	_sim->_spawnLocations.push_back(sf::Vector2f(SolidParticle::_size * 4, SolidParticle::_size * 6.3));
	_sim->_spawnLocations.push_back(sf::Vector2f(SolidParticle::_size * 4, SolidParticle::_size * 7.6));

	_sim->_spawnVelocities.push_back(sf::Vector2f(400, 0));
	_sim->_spawnVelocities.push_back(sf::Vector2f(400, 0));
	_sim->_spawnVelocities.push_back(sf::Vector2f(400, 0));
	
	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = false;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
	_sim->_spawnParticles = true;

	_sim->_spawnDelay = 0.005;
	_sim->_maxNumParticles = 10000;
}


// _________________________________________________________________________________
void InitManager::init_rain_simulation(int size, int zoom) {

	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 3);

	// Add boundaries to the side
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> line1 = placeParticleLine(pos, pos + sf::Vector2f(0, size * 2 + 6), solid, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}
	pos = sf::Vector2f(size * 2, 0);
	line1 = placeParticleLine(pos, pos + sf::Vector2f(0, size * 2 + 6), solid, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}

	
	// Add Box in the middle
	pos = sf::Vector2f(size - 10, SolidParticle::_size * 40);
	line1 = placeBox(pos, 10, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}
	pos = sf::Vector2f(size - 8, SolidParticle::_size * 41);
	line1 = placeBox(pos, 8, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}
	pos = sf::Vector2f(size - 6, SolidParticle::_size * 42);
	line1 = placeBox(pos, 6, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}
	pos = sf::Vector2f(size - 4, SolidParticle::_size * 43);
	line1 = placeBox(pos, 4, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}
	pos = sf::Vector2f(size - 2, SolidParticle::_size * 44);
	line1 = placeBox(pos, 2, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}

	pos = sf::Vector2f(0, 0);
	for (int i = 2; i < size - 1; i++) {
		_sim->_spawnLocations.push_back(sf::Vector2f(SolidParticle::_size * i, 0));
		_sim->_spawnVelocities.push_back(sf::Vector2f(0, 400));
	}

	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = false;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
	_sim->_spawnParticles = true;

	_sim->_spawnDelay = 0.0037;
	_sim->_maxNumParticles = 10000;
}

// _________________________________________________________________________________
void InitManager::init_rain2_simulation(int size, int zoom) {

	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 3);

	// Add boundaries to the side
	sf::Vector2f pos = sf::Vector2f(0, 0);
	std::vector<Particle> line1 = placeParticleLine(pos, pos + sf::Vector2f(0, size * 2 + 6), solid, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}
	pos = sf::Vector2f(size * 2, 0);
	line1 = placeParticleLine(pos, pos + sf::Vector2f(0, size * 2 + 6), solid, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}


	// Add first box on the left
	pos = sf::Vector2f(size - 22, SolidParticle::_size * 50);
	line1 = placeBox(pos, 10, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}
	pos = sf::Vector2f(size - 20, SolidParticle::_size * 51);
	line1 = placeBox(pos, 8, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}
	// Add second Box on the right
	pos = sf::Vector2f(size, SolidParticle::_size * 50);
	line1 = placeBox(pos, 10, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}
	pos = sf::Vector2f(size + 2, SolidParticle::_size * 51);
	line1 = placeBox(pos, 8, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}

	// Add third Box in the middle
	pos = sf::Vector2f(size - 22, SolidParticle::_size * 40);
	line1 = placeBox(pos, 10, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}
	pos = sf::Vector2f(size - 20, SolidParticle::_size * 41);
	line1 = placeBox(pos, 8, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}

	// Add triangle in the middle
	pos = sf::Vector2f(size - 22, SolidParticle::_size * 39);
	line1 = placeTriangle(pos, pos + sf::Vector2f(20, 0), pos + sf::Vector2f(10, -40), solid, _sim->_particles.size());
	for (int i = 0; i < line1.size(); i++) {
		_sim->_particles.push_back(line1[i]);
	}

	pos = sf::Vector2f(0, 0);
	for (int i = 2; i < size - 1; i++) {
		_sim->_spawnLocations.push_back(sf::Vector2f(SolidParticle::_size * i, 0));
		_sim->_spawnVelocities.push_back(sf::Vector2f(0, 400));
	}

	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = false;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
	_sim->_spawnParticles = true;

	_sim->_spawnDelay = 0.0037;
	_sim->_maxNumParticles = 10000;
}

// ______________________________________________________________________
void InitManager::init_fountain_simulation() {
	int size = 500;
	_sim->_zoomFactor = 1;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);

	std::vector<Particle> fountain = place_fountain(sf::Vector2f(50, 800), _sim);
	for (int i = 0; i < fountain.size(); i++) {
		_sim->_particles.push_back(fountain[i]);
	}

	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = false;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
	_sim->_spawnParticles = true;

	_sim->_maxNumParticles = 20000;
}

// _________________________________________________________________________________
void InitManager::init_moving_wall_simulation(int size, int zoom) {
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
	float xOffsetSize = 0.02;
	float yOffsetSize = 0.02;
	pos = sf::Vector2f(SolidParticle::_size * 2 + xOffsetSize / 2, SolidParticle::_size * (size - 2));
	for (int i = 0; i < size / 4; i++) {
		for (int j = 0; j < size - 3; j++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		xOffsetSize = xOffsetSize * -1;
		pos.x = SolidParticle::_size * 2 + xOffsetSize;
		pos.y -= FluidParticle::_size + yOffsetSize;
	}

	// Add moving wall
	MovingObject object1 = MovingObject();
	pos = sf::Vector2f(SolidParticle::_size * size / 2, SolidParticle::_size * size / 2);
	std::vector<Particle> wall = placeBox(pos, 10, solid);
	for (int i = 0; i < wall.size(); i++) {
		_sim->_particles.push_back(wall[i]);
		Particle* movingParticle = &_sim->_particles.at(_sim->_particles.size() - 1);
		movingParticle->_type = moving;
		movingParticle->_velocity = sf::Vector2f(0, 30);
		object1._particles.push_back(movingParticle);
		object1._ids.push_back(movingParticle->_id);
	}
	pos = pos + sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	wall = placeBox(pos, 8, solid);
	for (int i = 0; i < wall.size(); i++) {
		_sim->_particles.push_back(wall[i]);
		Particle* movingParticle = &_sim->_particles.at(_sim->_particles.size() - 1);
		movingParticle->_type = moving;
		movingParticle->_velocity = sf::Vector2f(0, 30);
		object1._particles.push_back(movingParticle);
		object1._ids.push_back(movingParticle->_id);
	}
	
	object1._conditions.push_back(sf::Vector2f(100000000, SolidParticle::_size * (size - 3)));
	object1._conditions.push_back(sf::Vector2f(0, SolidParticle::_size * size / 2));
	object1._conditionBigger.push_back(true);
	object1._conditionBigger.push_back(false);
	object1._directions.push_back(sf::Vector2f(0, 30));
	object1._directions.push_back(sf::Vector2f(0, -30));
	object1._state = 0;
	_sim->_movingObjects.push_back(object1);



	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = false;
	_sim->_spawnParticles = false;
	_sim->_moveSolids = true;
}

// _________________________________________________________________________________
void InitManager::init_wave_generator_simulation(int size, int zoom) {
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
	float xOffsetSize = 0.02;
	float yOffsetSize = 0.02;
	pos = sf::Vector2f(SolidParticle::_size * 2 + xOffsetSize / 2, SolidParticle::_size * (size - 2));
	for (int i = 0; i < size / 4; i++) {
		for (int j = 0; j < size - 3; j++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		xOffsetSize = xOffsetSize * -1;
		pos.x = SolidParticle::_size * 2 + xOffsetSize;
		pos.y -= FluidParticle::_size + yOffsetSize;
	}

	// Add moving wall
	MovingObject object1 = MovingObject();
	pos = sf::Vector2f(SolidParticle::_size * -size * 1 / 4, SolidParticle::_size * size * 3/4);
	std::vector<Particle> wall = placeBox(pos, SolidParticle::_size * size * 1 / 8, solid);
	for (int i = 0; i < wall.size(); i++) {
		_sim->_particles.push_back(wall[i]);
		Particle* movingParticle = &_sim->_particles.at(_sim->_particles.size() - 1);
		movingParticle->_type = moving;
		movingParticle->_velocity = sf::Vector2f(size, 0);
		object1._particles.push_back(movingParticle);
		object1._ids.push_back(movingParticle->_id);
	}
	pos = pos + sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	wall = placeBox(pos, SolidParticle::_size * (size * 1 / 8 - 1), solid);
	for (int i = 0; i < wall.size(); i++) {
		_sim->_particles.push_back(wall[i]);
		Particle* movingParticle = &_sim->_particles.at(_sim->_particles.size() - 1);
		movingParticle->_type = moving;
		movingParticle->_velocity = sf::Vector2f(size, 0);
		object1._particles.push_back(movingParticle);
		object1._ids.push_back(movingParticle->_id);
	}

	object1._conditions.push_back(sf::Vector2f(SolidParticle::_size * size * 1 / 4, size * 100000));
	object1._conditions.push_back(sf::Vector2f(SolidParticle::_size * -size * 1 / 4, 0));
	object1._conditionBigger.push_back(true);
	object1._conditionBigger.push_back(false);
	object1._directions.push_back(sf::Vector2f(size, 0));
	object1._directions.push_back(sf::Vector2f(-size, 0));
	object1._state = 0;
	_sim->_movingObjects.push_back(object1);



	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = false;
	_sim->_spawnParticles = false;
	_sim->_moveSolids = true;
}


// _________________________________________________________________________________
void  InitManager::init_cube_drop_simulation(int size, int zoom) {

	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);

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

	float xOffsetSize = 0.02;	// 0.3
	float yOffsetSize = 0.02;	// 0.05
	float xFluidSize = 1;
	float yFluidSize = 1;
	xFluidSize = size / 2; yFluidSize = size / 2;
	pos = sf::Vector2f(SolidParticle::_size * size / 2, SolidParticle::_size * size / 1.3);
	pos.x = SolidParticle::_size * size / 2 + xOffsetSize - xFluidSize;
	for (int i = 0; i < yFluidSize; i++) {
		for (int j = 0; j < xFluidSize; j++) {
			_sim->_particles.push_back(FluidParticle(_sim->_particles.size(), pos));
			pos.x += FluidParticle::_size;
		}
		xOffsetSize = xOffsetSize * -1;
		pos.x = SolidParticle::_size * size / 2 + xOffsetSize - xFluidSize;
		pos.y -= FluidParticle::_size + yOffsetSize;
	}
	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = true;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
	_sim->_spawnParticles = false;

}

// _____________________________________________________________________________________________________________________

void InitManager::init_complex_2_simulation(int size, int zoom) {
	_sim->_zoomFactor = zoom;
	_sim->_renderer = Renderer(_sim->_zoomFactor, FluidParticle::_size, SolidParticle::_size, _sim->_neighborRadius);
	_sim->_hashManager = CompactHashManager(_sim->_neighborRadius, size / 2, size / 2);

	sf::Vector2f lowerLeft = sf::Vector2f(0, size * SolidParticle::_size);
	sf::Vector2f pos1 = lowerLeft;
	sf::Vector2f pos2 = lowerLeft;
	sf::Vector2f pos3 = lowerLeft;
	sf::Vector2f pos4 = lowerLeft;
	std::vector<Particle> wall;

	// Place pond
	pos1.y -= size * SolidParticle::_size * 0.3;
	for (int i = 0; i < size * 0.3; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos1));
		pos1.y += SolidParticle::_size;
	}
	for (int i = 0; i < size; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos1));
		pos1.x += SolidParticle::_size;
	}
	for (int i = 0; i < size * 0.3; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos1));
		pos1.y -= SolidParticle::_size;
	}

	// Place water pool
	pos1.y = size * SolidParticle::_size * 0.1;
	pos1.x = size * SolidParticle::_size * 0.1;
	for (int i = 0; i < size * 0.3; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos1));
		pos1.y += SolidParticle::_size;
	}
	for (int i = 0; i < size * 0.1; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos1));
		pos1.x += SolidParticle::_size;
	}
	wall = placeParticleLine(pos1, pos1 + sf::Vector2f(SolidParticle::_size * size * 0.05, SolidParticle::_size * size * 0.03),
		solid, _sim->_particles.size());
	for (int i = 0; i < wall.size(); i++) {
		_sim->_particles.push_back(wall[i]);
	}
	pos1 += sf::Vector2f(SolidParticle::_size * size * 0.05, SolidParticle::_size * size * 0.03);

	pos2 = pos1;
	pos1.x += SolidParticle::_size * size * 0.05;
	pos3 = pos1;

	wall = placeParticleLine(pos1, pos1 + sf::Vector2f(SolidParticle::_size * size * 0.05, -SolidParticle::_size * size * 0.03),
		solid, _sim->_particles.size());
	for (int i = 0; i < wall.size(); i++) {
		_sim->_particles.push_back(wall[i]);
	}
	pos1 += sf::Vector2f(SolidParticle::_size * size * 0.05, -SolidParticle::_size * size * 0.03);

	for (int i = 0; i < size * 0.1; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos1));
		pos1.x += SolidParticle::_size;
	}
	for (int i = 0; i < size * 0.3; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos1));
		pos1.y -= SolidParticle::_size;
	}

	// Place tube
	pos1 = pos2;
	for (int i = 0; i < size * 0.15; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos1));
		pos1.y += SolidParticle::_size;
	}
	pos2 = pos1;
	pos1 = pos3;
	for (int i = 0; i < size * 0.15; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos1));
		pos1.y += SolidParticle::_size;
	}
	pos3 = pos1;

	// + curvature
	int radius = size * SolidParticle::_size * 0.1;
	float phi = 3.14;
	pos4 = pos3 + sf::Vector2f(radius, -std::sqrt(radius));
	pos4.y += SolidParticle::_size * 3.5;
	for (int i = 0; i < size * 0.3; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos4 + sf::Vector2f(radius * std::cos(phi), radius * std::sin(phi))));
		phi -= 0.039;
	}
	pos3 = pos4 + sf::Vector2f(radius * std::cos(phi), radius * std::sin(phi));
	radius = size * SolidParticle::_size * 0.1 + SolidParticle::_size * size * 0.05;
	phi = 3.14;
	for (int i = 0; i < size * 0.47; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos4 + sf::Vector2f(radius * std::cos(phi), radius * std::sin(phi))));
		phi -= 0.024;
	}
	pos2 = pos4 + sf::Vector2f(radius * std::cos(phi), radius * std::sin(phi));


	// Add moving wall
	MovingObject object1 = MovingObject();
	pos1 = pos2;
	wall = placeParticleLine(pos2, pos3,
		solid, _sim->_particles.size());
	for (int i = 0; i < wall.size(); i++) {
		_sim->_particles.push_back(wall[i]);
		Particle* movingParticle = &_sim->_particles.at(_sim->_particles.size() - 1);
		movingParticle->_type = moving;
		movingParticle->_velocity = sf::Vector2f((pos2.x - pos3.x) * 1, (pos2.y - pos3.y) * 1);
		object1._particles.push_back(movingParticle);
		object1._ids.push_back(movingParticle->_id);
	}
	
	object1._conditions.push_back(sf::Vector2f(1.07 * pos2.x, SolidParticle::_size * 20 * size));
	object1._conditions.push_back(sf::Vector2f(1.08 * pos2.x, SolidParticle::_size * 20 * size));
	object1._conditions.push_back(sf::Vector2f(0.93 * pos2.x, 0));
	object1._conditions.push_back(sf::Vector2f(0.92 * pos2.x, 0));
	object1._conditionBigger.push_back(true);
	object1._conditionBigger.push_back(true);
	object1._conditionBigger.push_back(false);
	object1._conditionBigger.push_back(false);
	object1._directions.push_back(sf::Vector2f((pos2.x - pos3.x) * 1, (pos2.y - pos3.y) * 1));
	object1._directions.push_back(sf::Vector2f((pos2.x - pos3.x) * 0.1, (pos2.y - pos3.y) * 0.1));
	object1._directions.push_back(sf::Vector2f(-(pos2.x - pos3.x) * 0.5, -(pos2.y - pos3.y) * 0.5));
	object1._directions.push_back(sf::Vector2f(-(pos2.x - pos3.x) * 0.05, -(pos2.y - pos3.y) * 0.05));
	object1._state = 0;
	_sim->_movingObjects.push_back(object1);





	// Place side spawns
	pos1 = sf::Vector2f(0, SolidParticle::_size * 4);
	for (int i = 0; i < size * 0.055; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos1));
		pos1.x += SolidParticle::_size;
	}
	pos1 = sf::Vector2f(0, SolidParticle::_size * 9.4);
	for (int i = 0; i < size * 0.055; i++) {
		_sim->_particles.push_back(SolidParticle(_sim->_particles.size(), pos1));
		pos1.x += SolidParticle::_size;
	}
	_sim->_spawnLocations.push_back(sf::Vector2f(SolidParticle::_size * 4, SolidParticle::_size * 5));
	_sim->_spawnLocations.push_back(sf::Vector2f(SolidParticle::_size * 4, SolidParticle::_size * 6.1));
	_sim->_spawnLocations.push_back(sf::Vector2f(SolidParticle::_size * 4, SolidParticle::_size * 7.2));
	_sim->_spawnLocations.push_back(sf::Vector2f(SolidParticle::_size * 4, SolidParticle::_size * 8.3));
	
	_sim->_spawnVelocities.push_back(sf::Vector2f(400, 0));
	_sim->_spawnVelocities.push_back(sf::Vector2f(400, 0));
	_sim->_spawnVelocities.push_back(sf::Vector2f(400, 0));
	_sim->_spawnVelocities.push_back(sf::Vector2f(400, 0));
	
	_sim->_spawnDelay = 0.005;
	_sim->_maxNumParticles = 10000;

	// Add moving cube
	MovingObject object2 = MovingObject();
	pos1 = sf::Vector2f(SolidParticle::_size * size * 0.75, SolidParticle::_size * size * 0.8);
	std::vector<Particle> box = placeBox(pos1, 40, _sim->_particles.size());
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
		Particle* movingParticle = &_sim->_particles.at(_sim->_particles.size() - 1);
		movingParticle->_type = moving;
		movingParticle->_velocity = sf::Vector2f(0, size / 4);
		object2._particles.push_back(movingParticle);
		object2._ids.push_back(movingParticle->_id);
	}
	pos1 = pos1 + sf::Vector2f(SolidParticle::_size, SolidParticle::_size);
	box = placeBox(pos1, 38, _sim->_particles.size());
	for (int i = 0; i < box.size(); i++) {
		_sim->_particles.push_back(box[i]);
		Particle* movingParticle = &_sim->_particles.at(_sim->_particles.size() - 1);
		movingParticle->_type = moving;
		movingParticle->_velocity = sf::Vector2f(0, size / 4);
		object2._particles.push_back(movingParticle);
		object2._ids.push_back(movingParticle->_id);
	}
	
	object2._conditions.push_back(sf::Vector2f(SolidParticle::_size * 20 * size, lowerLeft.y - 5 * SolidParticle::_size));
	object2._conditions.push_back(sf::Vector2f(0, SolidParticle::_size * size * 0.5));
	std::cout << object2._conditions[0].x << " " << object2._conditions[0].y << "	 " << object2._conditions[1].x << " " << object2._conditions[1].y << std::endl;
	object2._conditionBigger.push_back(true);
	object2._conditionBigger.push_back(false);
	object2._directions.push_back(sf::Vector2f(0, size / 4));
	object2._directions.push_back(sf::Vector2f(0, -size / 4));
	object2._state = 0;
	_sim->_movingObjects.push_back(object2);

	


	_sim->_moveParticles = true;
	_sim->_testNeighbors = false;
	_sim->_testKernel = false;
	_sim->_printFPS = false;
	_sim->_printParticleInfo = false;
	_sim->_deleteParticles = true;
	_sim->_spawnParticles = true;
	_sim->_moveSolids = true;
}