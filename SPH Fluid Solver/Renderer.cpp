// Paul Kull, 2022

#define _USE_MATH_DEFINES

#include <iostream>

#include "./Renderer.h"
#include "./Parameters.h"


static sf::Font font;

// ___________________________________________________________
Renderer::Renderer(float zoomFactor, float fluidSize, float solidSize, float searchRadius) {
	_zoomFactor = zoomFactor;

	_fluidParticleShapes = std::vector<sf::CircleShape>();
	_solidParticleShapes = std::vector<sf::CircleShape>();
	_watchedParticleShape = sf::CircleShape();
	_arrowBodies = std::vector<sf::RectangleShape>();
	_arrowHeads = std::vector<sf::CircleShape>();
	_infoPanel = sf::RectangleShape(sf::Vector2f(300, 800));
	_graphBackground = sf::RectangleShape(sf::Vector2f(250, 100));
	_graphShapes = std::vector<sf::CircleShape>();
	_fluidParticleShapes.clear();
	_solidParticleShapes.clear();

	_fluidShapeRadius = _zoomFactor * fluidSize / 2;
	_solidShapeRadius = _zoomFactor * fluidSize / 2;

	_watchedParticleShape.setFillColor(sf::Color::Green);
	_watchedParticleShape.setRadius(_fluidShapeRadius);
	_watchedParticleShape.setPosition(sf::Vector2f(-100, -100));

	float outlineThickness = 3;
	_searchRadiusShape = sf::CircleShape();
	_searchRadiusShape.setRadius(searchRadius * _zoomFactor);
	_searchRadiusShape.setFillColor(sf::Color::Transparent);
	_searchRadiusShape.setOutlineColor(sf::Color::Magenta);
	_searchRadiusShape.setOutlineThickness(outlineThickness);
	_searchRadiusShape.setPosition(sf::Vector2f(-100 * _zoomFactor, -100 * _zoomFactor));
	_searchRadiusOffset.x = (- searchRadius + FluidParticle::_size / 2) * _zoomFactor;
	_searchRadiusOffset.y = (- searchRadius + FluidParticle::_size / 2) * _zoomFactor;

	_infoPanel.setPosition(sf::Vector2f(800, 0));
	_infoPanel.setFillColor(sf::Color::Color(200, 200, 200));
	_graphBackground.setPosition(sf::Vector2f(825, 650));
	_graphBackground.setFillColor(sf::Color::Color(250, 250, 250));
	_graphBackground.setOutlineColor(sf::Color::Black);
	_graphBackground.setOutlineThickness(1);
	_graphShapes.clear();
	_graphShapesFull = false;


	if (!font.loadFromFile("../Resources/times new roman.ttf")) {
		throw("ERROR: Could not load font");
	}

	_description = sf::Text();
	_description.setFont(font);
	_information = sf::Text();
	_information.setFont(font);
	
	_description.setPosition(sf::Vector2f(820, 10));
	_description.setString(
		"Current Time:\n\nNumber of particles:\n\n# moving Particles\n\nUpdates per second:\n\nAverage Fluid Density : \n\nMaximum Timestep : \n\n\n\nCurrent Particle:\n\nDensity: ");
	_description.setCharacterSize(15);
	_description.setFillColor(sf::Color::Black);
	_description.setStyle(sf::Text::Underlined);

	_information.setCharacterSize(15);
	_information.setFillColor(sf::Color::Black);
	_information.setPosition(sf::Vector2f(1000, 10));
}


// ___________________________________________________________
void Renderer::update_graphics(std::vector<Particle>* particles, int numFluids, int watchedParticleId, std::vector<int> markedParticlesId, std::vector<int> testedParticlesId) {


	// Check whether or not we have the correct amount of shapes
	if (_fluidParticleShapes.size() + _solidParticleShapes.size() != particles->size()) {
		// Add shapes until there are enough
		while (_fluidParticleShapes.size() + _solidParticleShapes.size() < particles->size()) {
			if (particles->at(_fluidParticleShapes.size() + _solidParticleShapes.size())._type == solid) {
				_solidParticleShapes.push_back(sf::CircleShape(_solidShapeRadius));
				_solidParticleShapes.back().setPosition(particles->at(_fluidParticleShapes.size() + _solidParticleShapes.size() - 1)._position * _zoomFactor);
				_solidParticleShapes.back().setFillColor(sf::Color::White);
			}
			else if (particles->at(_fluidParticleShapes.size() + _solidParticleShapes.size())._type == fluid) {
				_fluidParticleShapes.push_back(sf::CircleShape(_fluidShapeRadius));
			}
		}
		// Delete shapes until there are little enough
		while (_fluidParticleShapes.size() + _solidParticleShapes.size() > particles->size()) {
			_fluidParticleShapes.pop_back();
			
		}
	}

	int numParticles = particles->size();
	int fluidIndex = 0;
	int speed = 0;

	// Update each fluid shapes position
	for (int i = 0; i < numParticles; i++) {
		if (particles->at(i)._id == watchedParticleId) {
			_watchedParticleShape.setPosition(particles->at(i)._position * _zoomFactor);
			_searchRadiusShape.setPosition(particles->at(i)._position * _zoomFactor + _searchRadiusOffset);
			update_arrows(&particles->at(i));
		}
		if (particles->at(i)._type == solid) { continue; }

		_fluidParticleShapes[fluidIndex].setPosition(particles->at(i)._position * _zoomFactor);


		int density = std::min((int)(particles->at(i)._density * 100), 255);
		speed = std::min((int)(Functions::calculate_distance_norm(particles->at(i)._velocity) * 0.6), 255);
		_fluidParticleShapes[fluidIndex].setFillColor(FluidParticle::_stasisColor + sf::Color::Color(density, speed, 0));

		int numTestedParticles = testedParticlesId.size();
		for (int j = 0; j < numTestedParticles; j++) {
			if (particles->at(i)._id == testedParticlesId[j]) {
				_fluidParticleShapes[fluidIndex].setFillColor(sf::Color::Cyan);
			}
		}

		int numMarkedParticles = markedParticlesId.size();
		for (int j = 0; j < numMarkedParticles; j++) {
			if (particles->at(i)._id == markedParticlesId[j]) {
				_fluidParticleShapes[fluidIndex].setFillColor(sf::Color::Red);
			}
		}
		fluidIndex++;
	}
}


// ___________________________________________________________
void Renderer::update_information(sf::Time time, int numParticles, int numFluidParticles, float numUpdates, float avgDensity, float maxStep, float watchedParticleDensity) {
	// Information in the box
	_timeInfo = std::to_string(time.asSeconds());
	_numParticlesInfo = std::to_string(numParticles);
	_numFluidsInfo = std::to_string(numFluidParticles);
	_numUpdatesInfo = std::to_string(numUpdates);
	_avgDensityInfo = std::to_string(avgDensity);
	_maxStepInfo = std::to_string(maxStep);
	_watchedParticleDensity = std::to_string(watchedParticleDensity);

	_timeInfo.resize(4, ' ');
	_numUpdatesInfo.resize(3, ' ');
	_avgDensityInfo.resize(6, ' ');
	_maxStepInfo.resize(6, ' ');
	_watchedParticleDensity.resize(4, ' ');

	_information.setString(_timeInfo + "\n\n" + _numParticlesInfo + "\n\n" + _numFluidsInfo + "\n\n" + _numUpdatesInfo + "\n\n"
		+ _avgDensityInfo + "\n\n" + _maxStepInfo + "\n\n\n\n\n\n" + _watchedParticleDensity);


	// Take care of the Graph
	for (int i = 0; i < _graphShapes.size(); i++) {
		_graphShapes[i].move(-Parameters::GRAPH_SPEED, 0);
		if (_graphShapes[i].getPosition().x < _graphBackground.getPosition().x) {
			_graphShapes[i].setPosition(sf::Vector2f(1000, 750 - 50 * std::pow(avgDensity, Parameters::GRAPH_ZOOM)));
			_graphShapesFull = true;
		}
	}
	if (!_graphShapesFull) {
		_graphShapes.push_back(sf::CircleShape(2));
		_graphShapes.back().setFillColor(sf::Color::Blue);
		_graphShapes.back().setPosition(sf::Vector2f(1000, 750 - 50 * std::pow(avgDensity, Parameters::GRAPH_ZOOM)));
	}
}

// ___________________________________________________________
void Renderer::update_arrows(Particle* watchedParticle) {
	float scalingFactor = 0.05;
	int thickness = 5;

	_arrowBodies.clear();

	float sizeAccel = Functions::calculate_distance_norm(watchedParticle->_acceleration) * scalingFactor;
	float sizePress = Functions::calculate_distance_norm(watchedParticle->_pressureAcc) * scalingFactor;

	// Calculate angle difference to vector (1, 0)
	float angleAccel = std::acos(watchedParticle->_acceleration.x /
		Functions::calculate_distance_norm(watchedParticle->_acceleration)) * 180 / M_PI;
	float anglePress = std::acos(watchedParticle->_pressureAcc.x /
		Functions::calculate_distance_norm(watchedParticle->_pressureAcc)) * 180 / M_PI;;

	sf::Vector2f particleSizeOffsetAccel = sf::Vector2f(FluidParticle::_size / 2 + thickness / _zoomFactor / 2, FluidParticle::_size / 2);
	sf::Vector2f particleSizeOffsetPress = sf::Vector2f(FluidParticle::_size / 2 + thickness / _zoomFactor / 2, FluidParticle::_size / 2);
	// Add arrow for acceleration
	_arrowBodies.push_back(sf::RectangleShape(sf::Vector2f(sizeAccel, thickness)));
	_arrowBodies.back().setFillColor(sf::Color::Yellow);
	_arrowBodies.back().setPosition((watchedParticle->_position + particleSizeOffsetAccel) * _zoomFactor);
	if (watchedParticle->_acceleration.y > 0) {
		_arrowBodies.back().rotate(angleAccel);
	}
	else {
		_arrowBodies.back().rotate(-angleAccel);
	}
	// Add arrow for pressure
	_arrowBodies.push_back(sf::RectangleShape(sf::Vector2f(sizePress, thickness)));
	_arrowBodies.back().setFillColor(sf::Color::Color(150, 200, 100));
	_arrowBodies.back().setPosition((watchedParticle->_position + particleSizeOffsetPress) * _zoomFactor);
	if (watchedParticle->_pressureAcc.y > 0) {
		_arrowBodies.back().rotate(anglePress);
	}
	else {
		_arrowBodies.back().rotate(-anglePress);
	}

}



// ___________________________________________________________
void Renderer::draw(sf::RenderWindow* window) {

	window->clear(sf::Color::Black);

	for (int i = 0; i < _solidParticleShapes.size(); i++) {
		window->draw(_solidParticleShapes[i]);
	}
	for (int i = 0; i < _fluidParticleShapes.size(); i++) {
		window->draw(_fluidParticleShapes[i]);
	}
	for (int i = 0; i < _arrowBodies.size(); i++) {
		window->draw(_arrowBodies[i]);
	}
	window->draw(_watchedParticleShape);
	// window->draw(_searchRadiusShape);
	window->draw(_infoPanel);
	window->draw(_graphBackground);
	for (int i = 0; i < _graphShapes.size(); i++) {
		window->draw(_graphShapes[i]);
	}
	window->draw(_description);
	window->draw(_information);

	window->display();
}
