// Paul Kull, 2022

#include "./Renderer.h"


// ___________________________________________________________
Renderer::Renderer(int fluidRadius, int solidRadius) {
	_particleShapes = std::vector<sf::CircleShape>();
	_particleShapes.clear();

	_fluidParticleRadius = fluidRadius;
	_solidParticleRadius = fluidRadius;
}


// ___________________________________________________________
void Renderer::update_graphics(std::vector<Particle>* particles) {


	// Check whether or not we have the correct amount of shapes
	if (_particleShapes.size() != particles->size()) {
		// Add shapes until there are enough
		while (_particleShapes.size() < particles->size()) {
			_particleShapes.push_back(sf::CircleShape());
		}
		// Delete shapes until there are little enough
		while (_particleShapes.size() > particles->size()) {
			_particleShapes.pop_back();
		}
	}

	int numParticles = _particleShapes.size();

	// Update each shapes position
	for (int i = 0; i < numParticles; i++) {
		_particleShapes[i].setPosition(particles->at(i)._position);
		switch (particles->at(i)._type) {
		case(fluid):
			_particleShapes[i].setRadius(_fluidParticleRadius);
			break;
		case(solid):
			_particleShapes[i].setRadius(_solidParticleRadius);
			break;
		}
		_particleShapes[i].setFillColor(particles->at(i)._stasisColor);
	}
}



// ___________________________________________________________
void Renderer::draw(sf::RenderWindow* window) {

	window->clear(sf::Color::Black);
	int numShapes = _particleShapes.size();

	for (int i = 0; i < numShapes; i++) {
		window->draw(_particleShapes[i]);
	}

	window->display();
}