// Paul Kull, 2022

#include <iostream>

#include "./Renderer.h"


// ___________________________________________________________
Renderer::Renderer(float zoomFactor, int fluidRadius, int solidRadius, int searchRadius) {
	_zoomFactor = zoomFactor;

	_particleShapes = std::vector<sf::CircleShape>();
	_particleShapes.clear();

	int outlineThickness = 1;
	_searchRadiusShape = sf::CircleShape();
	_searchRadiusShape.setRadius(searchRadius * _zoomFactor);
	_searchRadiusShape.setFillColor(sf::Color::Transparent);
	_searchRadiusShape.setOutlineColor(sf::Color::Magenta);
	_searchRadiusShape.setOutlineThickness(outlineThickness);
	_searchRadiusShape.setPosition(sf::Vector2f(-100 * _zoomFactor, -100 * _zoomFactor));
	_searchRadiusOffset.x = ( - searchRadius + outlineThickness)* _zoomFactor;
	_searchRadiusOffset.y = ( - searchRadius + outlineThickness) * _zoomFactor;

	_fluidShapeRadius = _zoomFactor * fluidRadius;
	_solidShapeRadius = _zoomFactor * fluidRadius;

	if (!_font.loadFromFile("../Resources/times new roman.ttf")) {
		std::cout << "ERROR: Could not load font";
	}
	// _description = sf::Text();
	// _description.setFont(_font);
	// _information = sf::Text();
	// _information.setFont(_font);
	// 
	// _description.setString("# of particles:\nUpdates per second:");
	// _description.setCharacterSize(10);
	// _description.setFillColor(sf::Color::White);
	// // _description.setStyle(sf::Text::Underlined);
	// _information.setPosition(sf::Vector2f(10, 10));
	// _information.setCharacterSize(10);
	// _information.setFillColor(sf::Color::White);
	// _information.setPosition(sf::Vector2f(30, 10));
}


// ___________________________________________________________
void Renderer::update_graphics(std::vector<Particle>* particles, int watchedParticleId, std::vector<int> markedParticlesId) {


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
		_particleShapes[i].setPosition(particles->at(i)._position * _zoomFactor);
		switch (particles->at(i)._type) {
		case(fluid):
			_particleShapes[i].setRadius(_fluidShapeRadius);
			break;
		case(solid):
			_particleShapes[i].setRadius(_solidShapeRadius);
			break;
		}

		_particleShapes[i].setFillColor(particles->at(i)._stasisColor);

		if (particles->at(i)._id == watchedParticleId) {
			_particleShapes[i].setFillColor(sf::Color::Green);
			_searchRadiusShape.setPosition(particles->at(i)._position * _zoomFactor + _searchRadiusOffset);
			continue;
		}
		int numMarkedParticles = markedParticlesId.size();
		for (int j = 0; j < numMarkedParticles; j++) {
			if (particles->at(i)._id == markedParticlesId[j]) {
				_particleShapes[i].setFillColor(sf::Color::Red);
			}
		}
	}
}


// ___________________________________________________________
void Renderer::update_information(int numParticles, int numUpdates) {
	// _information.setString(std::to_string(numParticles) + "\n" + std::to_string(numUpdates));
	std::cout << "numParticles: " << numParticles << "	numUpdates: " << numUpdates << "\n\n";
}

// ___________________________________________________________
void Renderer::draw(sf::RenderWindow* window) {

	window->clear(sf::Color::Black);
	int numShapes = _particleShapes.size();

	for (int i = 0; i < numShapes; i++) {
		window->draw(_particleShapes[i]);
	}
	window->draw(_searchRadiusShape);
	window->display();
}
