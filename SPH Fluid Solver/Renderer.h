#pragma once
// Paul Kull, 2022

#include <SFML/Graphics.hpp>
#include "./SolidParticle.h"
#include "./FluidParticle.h"
#include "./Functions.h"

class Renderer {
  public:
	std::vector<sf::CircleShape> _fluidParticleShapes;
	std::vector<sf::CircleShape> _solidParticleShapes;
	sf::CircleShape _watchedParticleShape;
	std::vector<sf::RectangleShape>_arrowBodies;
	std::vector<sf::CircleShape> _arrowHeads;
	sf::RectangleShape _infoPanel;
	sf::RectangleShape _graphBackground;
	std::vector<sf::CircleShape> _graphShapes;
	sf::CircleShape _searchRadiusShape;
	sf::Vector2f _searchRadiusOffset;
	sf::Text _description;
	sf::Text _information;

	std::string _timeInfo;
	std::string _numParticlesInfo;
	std::string _numFluidsInfo;
	std::string _numUpdatesInfo;
	std::string _avgDensityInfo;
	std::string _maxStepInfo;
	std::string _watchedParticleDensity;


	bool _graphShapesFull;

	float _zoomFactor;
	int _fluidShapeRadius;
	int _solidShapeRadius;

	Renderer(float zoomFactor = 1, float fluidSize = 1, float solidSize = 1, float searchRadius = 3);
	void update_graphics(std::vector<Particle>* particles, int numFluids, int watchedParticleId, std::vector<int> markedParticlesId, std::vector<int> testedParticlesId);
	void update_information(float time, int numParticles, int numFluidParticles, float numUpdates, float avgDensity = -1, float maxStep = -1, float watchedParticleDensity = -1);
	void update_arrows(Particle* watchedPartile);
	void draw(sf::RenderWindow* window);
};