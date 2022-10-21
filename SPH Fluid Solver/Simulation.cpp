// Paul Kull, 2022

#include "./Simulation.h"

// _________________________________________________________________________________
Simulation::Simulation(SimulationPreset preset = StuffedBox) {

	_particles = std::vector<Particle>();
	_particles.clear();

	switch (preset) {
	case Empty:
		init_empty_simulation();
		break;
	case StuffedBox:
		init_stuffed_box_simulation();
	}

}

// _________________________________________________________________________________
void Simulation::init_empty_simulation() {
	sf::VideoMode _videoMode = sf::VideoMode();
	_videoMode.size = sf::Vector2u(200, 200);
	_window.create(_videoMode, "SPH Fluid Solver");
	renderer = Renderer();
}

// _________________________________________________________________________________
void  Simulation::init_stuffed_box_simulation() {

	sf::VideoMode _videoMode = sf::VideoMode();
	_videoMode.size = sf::Vector2u(800, 800);
	_window.create(_videoMode, "SPH Fluid Solver");
	renderer = Renderer(5, 5);

	// Add Particles
	sf::Vector2f pos = sf::Vector2f(150, 200);
	for (int i = 0; i < 50; i++) {
		_particles.push_back(SolidParticle(pos));
		pos.y += 10;
	}
	for (int i = 0; i < 50; i++) {
		_particles.push_back(SolidParticle(pos));
		pos.x += 10;
	}
	for (int i = 0; i < 50; i++) {
		_particles.push_back(SolidParticle(pos));
		pos.y -= 10;
	}
	for (int i = 0; i < 50; i++) {
		_particles.push_back(SolidParticle(pos));
		pos.x -= 10;
	}

	pos.x = 160;
	pos.y = 210;

	for (int i = 0; i < 49; i++) {
		for (int ii = 0; ii < 49; ii++) {
			_particles.push_back(FluidParticle(pos));
			pos.y += 10;
		}
		pos.y = 210;
		pos.x += 10;
	}
}


// _________________________________________________________________________________
void update_physics() {
	return;
}




// _________________________________________________________________________________
void Simulation::run() {
	bool runSimulation = true;

	while (runSimulation) {
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape)) {
			_window.close();
			runSimulation = false;
			break;
		}
		// update_physics();
		renderer.update_graphics(&_particles);
		renderer.draw(&_window);
	}
}