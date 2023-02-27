// Paul Kull, 2022

#define _USE_MATH_DEFINES

#include <iostream>

#include "./Renderer.h"
#include "./Parameters.h"


static sf::Font font;

// ___________________________________________________________
Renderer::Renderer(float zoomFactor, float fluidSize, float solidSize, float searchRadius) {
	_zoomFactor = zoomFactor;

	_solidShapes = std::vector<sf::CircleShape>();
	_fluidShape = sf::CircleShape();
	_arrowBodies = std::vector<sf::RectangleShape>();
	_arrowHeads = std::vector<sf::CircleShape>();
	_infoPanel = sf::RectangleShape(sf::Vector2f(300, Parameters::WINDOW_HEIGHT));
	_graphBackground = sf::RectangleShape(sf::Vector2f(250, 100));
	_graphShapes = std::vector<sf::CircleShape>();

	_fluidShapeRadius = _zoomFactor * fluidSize / 2;
	_solidShapeRadius = _zoomFactor * fluidSize / 2;

	_fluidShape.setRadius(_fluidShapeRadius);

	float outlineThickness = 3;
	_searchRadiusShape = sf::CircleShape();
	_searchRadiusShape.setRadius(searchRadius * _zoomFactor);
	_searchRadiusShape.setFillColor(sf::Color::Transparent);
	_searchRadiusShape.setOutlineColor(sf::Color::Magenta);
	_searchRadiusShape.setOutlineThickness(outlineThickness);
	_searchRadiusShape.setPosition(sf::Vector2f(-100 * _zoomFactor, -100 * _zoomFactor));
	_searchRadiusOffset.x = (- searchRadius + FluidParticle::_size / 2) * _zoomFactor;
	_searchRadiusOffset.y = (- searchRadius + FluidParticle::_size / 2) * _zoomFactor;

	_infoPanel.setPosition(sf::Vector2f(Parameters::WINDOW_WIDTH - 300, 0));
	_infoPanel.setFillColor(sf::Color::Color(200, 200, 200));
	_graphBackground.setPosition(sf::Vector2f(Parameters::WINDOW_WIDTH - 275, 650));
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
	
	_description.setPosition(sf::Vector2f(Parameters::WINDOW_WIDTH - 280, 10));
	_description.setString(
		"Current Time:\n\nNumber of particles:\n\n# moving Particles\n\nUpdates per second:\n\nAverage Fluid Density : \n\nMaximum Velocity : \n\n\n\nCurrent Particle:\n\nDensity: ");
	_description.setCharacterSize(15);
	_description.setFillColor(sf::Color::Black);
	_description.setStyle(sf::Text::Underlined);

	_information.setCharacterSize(15);
	_information.setFillColor(sf::Color::Black);
	_information.setPosition(sf::Vector2f(Parameters::WINDOW_WIDTH - 100, 10));
}

// __________________________________________________________________________________
void Renderer::init_solids(std::vector<Particle>* particles) {
	for (int i = 0; i < particles->size(); i++) {
		if (particles->at(i)._type == solid) {
			sf::CircleShape newShape = sf::CircleShape();
			newShape.setFillColor(sf::Color::Color::White);
			newShape.setRadius(_solidShapeRadius);
			newShape.setPosition(particles->at(i)._position.x * _zoomFactor, particles->at(i)._position.y * _zoomFactor);
			_solidShapes.push_back(newShape);
		}
	}
}


// ___________________________________________________________
void Renderer::update_information(float time, int numParticles, int numFluidParticles, float numUpdates, float avgDensity,
	float maxVel, float watchedParticleDensity, bool updateGraph) {
	// Information in the box
	_timeInfo = std::to_string(time);
	_numParticlesInfo = std::to_string(numParticles);
	_numFluidsInfo = std::to_string(numFluidParticles);
	_numUpdatesInfo = std::to_string(numUpdates);
	_avgDensityInfo = std::to_string(avgDensity);
	_maxStepInfo = std::to_string(maxVel);
	_watchedParticleDensity = std::to_string(watchedParticleDensity);

	_timeInfo.resize(4, ' ');
	_numUpdatesInfo.resize(3, ' ');
	_avgDensityInfo.resize(6, ' ');
	_maxStepInfo.resize(6, ' ');
	_watchedParticleDensity.resize(4, ' ');

	_information.setString(_timeInfo + "\n\n" + _numParticlesInfo + "\n\n" + _numFluidsInfo + "\n\n" + _numUpdatesInfo + "\n\n"
		+ _avgDensityInfo + "\n\n" + _maxStepInfo + "\n\n\n\n\n\n" + _watchedParticleDensity);

	if (updateGraph) {
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
void Renderer::draw(sf::RenderWindow* window, std::vector<Particle>* particles,
	int watchedParticleId, std::vector<int> markedParticlesId,
	std::vector<int> testedParticlesId, bool updateArrows, bool drawGraph, bool drawArrows) {

	window->clear(sf::Color::Black);

	int numParticles = particles->size();

	// Update each fluid shapes position
	for (int i = 0; i < numParticles; i++) {
		if (particles->at(i)._type == solid) { continue; }

		_fluidShape.setPosition(particles->at(i)._position * _zoomFactor);
		_fluidShape.setFillColor(particles->at(i)._stasisColor + sf::Color::Color(0, particles->at(i)._colorFactor, 0));
		//if (particles->at(i)._type == fluid) {
		//	_particleShape.setRadius(_fluidShapeRadius);
		//}
		//else if (particles->at(i)._type == solid) {
		//	_particleShape.setRadius(_solidShapeRadius);
		//}
		//if (particles->at(i)._id == watchedParticleId) {
		//	_particleShape.setFillColor(sf::Color::Color::Red);
		//	if (updateArrows) {
		//		update_arrows(&particles->at(i));
		//	}
		//}
		window->draw(_fluidShape);
	}

	for (int i = 0; i < _solidShapes.size(); i++) {
		window->draw(_solidShapes[i]);
	}

	if (drawArrows) {
		for (int i = 0; i < _arrowBodies.size(); i++) {
			window->draw(_arrowBodies[i]);
		}
	}
	window->draw(_searchRadiusShape);
	window->draw(_infoPanel);
	if (drawGraph) {
		window->draw(_graphBackground);
		for (int i = 0; i < _graphShapes.size(); i++) {
			window->draw(_graphShapes[i]);
		}
	}
	window->draw(_description);
	window->draw(_information);

	window->display();
}
