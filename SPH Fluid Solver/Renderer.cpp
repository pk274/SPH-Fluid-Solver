// Paul Kull, 2022

#define _USE_MATH_DEFINES

#include <iostream>

#include "./Renderer.h"


// ___________________________________________________________
Renderer::Renderer(float zoomFactor, float fluidSize, float solidSize, float searchRadius) {
	_zoomFactor = zoomFactor;

	_particleShapes = std::vector<sf::CircleShape>();
	_arrowBodies = std::vector<sf::RectangleShape>();
	_arrowHeads = std::vector<sf::CircleShape>();
	_particleShapes.clear();

	_fluidShapeRadius = _zoomFactor * fluidSize / 2;
	_solidShapeRadius = _zoomFactor * fluidSize / 2;

	float outlineThickness = 2;
	_searchRadiusShape = sf::CircleShape();
	_searchRadiusShape.setRadius(searchRadius * _zoomFactor);
	_searchRadiusShape.setFillColor(sf::Color::Transparent);
	_searchRadiusShape.setOutlineColor(sf::Color::Magenta);
	_searchRadiusShape.setOutlineThickness(outlineThickness);
	_searchRadiusShape.setPosition(sf::Vector2f(-100 * _zoomFactor, -100 * _zoomFactor));
	_searchRadiusOffset.x = (- searchRadius + FluidParticle::_size / 2) * _zoomFactor;
	_searchRadiusOffset.y = (- searchRadius + FluidParticle::_size / 2) * _zoomFactor;

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
void Renderer::update_graphics(std::vector<Particle>* particles, int watchedParticleId, std::vector<int> markedParticlesId, std::vector<int> testedParticlesId) {


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

		int numTestedParticles = testedParticlesId.size();
		for (int j = 0; j < numTestedParticles; j++) {
			if (particles->at(i)._id == testedParticlesId[j]) {
				_particleShapes[i].setFillColor(sf::Color::Cyan);
			}
		}

		if (particles->at(i)._id == watchedParticleId) {
			_particleShapes[i].setFillColor(sf::Color::Green);
			_searchRadiusShape.setPosition(particles->at(i)._position * _zoomFactor + _searchRadiusOffset);
			std::cout << particles->at(i)._position.x << " " << particles->at(i)._position.y << "			"
				<< particles->at(i)._velocity.x << " " << particles->at(i)._velocity.y << "			" <<
				particles->at(i)._acceleration.x << " " << particles->at(i)._acceleration.y << std::endl;
			std::cout << "d: " << particles->at(i)._density << "			" << "p: " <<
				particles->at(i)._pressure << std::endl;
			update_arrows(particles, &particles->at(i));
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
void Renderer::update_information(int numParticles, float numUpdates) {
	// _information.setString(std::to_string(numParticles) + "\n" + std::to_string(numUpdates));
	std::cout << "\r                                   ";
	std::cout << "\rUpdates per Second: " << numUpdates << "	Number of Particles : " << numParticles << std::flush;
}

// ___________________________________________________________
void Renderer::update_arrows(std::vector<Particle>* particles, Particle* watchedParticle) {
	_arrowBodies.clear();
	_arrowHeads.clear();
	int scalingFactor = 10;
	float angleAccel = std::acos(watchedParticle->_acceleration.x /
		Functions::calculate_distance_norm(watchedParticle->_acceleration)) * 180 / M_PI;
	float sizeAccel = Functions::calculate_distance_norm(watchedParticle->_acceleration) * scalingFactor;
	float anglePress = std::acos(watchedParticle->_acceleration.x /
		Functions::calculate_distance_norm(watchedParticle->_pressureAcc)) * 180 / M_PI;;
	float sizePress = Functions::calculate_distance_norm(watchedParticle->_pressureAcc) * scalingFactor;
	int thickness = 10;
	sf::Vector2f particleSizeOffsetAccel = sf::Vector2f(FluidParticle::_size / 2 + thickness / _zoomFactor / 2, FluidParticle::_size / 2);
	sf::Vector2f particleSizeOffsetPress = sf::Vector2f(FluidParticle::_size / 2 - thickness / _zoomFactor / 2, FluidParticle::_size / 2);
	// Add arrow for acceleration
	_arrowBodies.push_back(sf::RectangleShape(sf::Vector2f(sizeAccel, thickness)));
	_arrowBodies[0].setFillColor(sf::Color::Yellow);
	_arrowBodies[0].setPosition((watchedParticle->_position + particleSizeOffsetAccel) * _zoomFactor);
	_arrowBodies[0].rotate(angleAccel);
	// Add arrow for pressure
	_arrowBodies.push_back(sf::RectangleShape(sf::Vector2f(sizePress, thickness)));
	_arrowBodies[1].setFillColor(sf::Color::Color(150, 200, 100));
	_arrowBodies[1].setPosition((watchedParticle->_position + particleSizeOffsetPress) * _zoomFactor);
	_arrowBodies[1].rotate(-anglePress);
	// _arrowHeads.push_back(sf::CircleShape(size / 10, 3));
	// _arrowHeads[0].setPosition((watchedParticle->_position + watchedParticle->_acceleration +
	// 	sf::Vector2f(FluidParticle::_size + thickness / 4 / _zoomFactor, - FluidParticle::_size)) * _zoomFactor);
	// _arrowHeads[0].setFillColor(sf::Color::Yellow);
	// _arrowHeads[0].rotate(angle + 90);


}


// ___________________________________________________________
void Renderer::draw(sf::RenderWindow* window) {

	window->clear(sf::Color::Black);
	int numShapes = _particleShapes.size();

	for (int i = 0; i < numShapes; i++) {
		window->draw(_particleShapes[i]);
	}
	for (int i = 0; i < _arrowHeads.size(); i++) {
		window->draw(_arrowHeads[i]);
	}
	for (int i = 0; i < _arrowBodies.size(); i++) {
		window->draw(_arrowBodies[i]);
	}
	window->draw(_searchRadiusShape);
	window->display();
}
