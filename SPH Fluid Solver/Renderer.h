#pragma once
// Paul Kull, 2022

#include <SFML/Graphics.hpp>
#include "./SolidParticle.h"
#include "./FluidParticle.h"
#include "./Functions.h"

class Renderer {
  public:
	std::vector<sf::CircleShape> _solidShapes;
	sf::CircleShape _fluidShape;
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

	std::string _applicationTimeInfo;
	std::string _simulatedTimeInfo;
	std::string _numParticlesInfo;
	std::string _numFluidsInfo;
	std::string _numUpdatesInfo;
	std::string _avgDensityInfo;
	std::string _maxStepInfo;
	std::string _avgSolverIters;
	std::string _watchedParticleDensity;

	sf::Texture _screenTexture;
	int _frameCounter;


	bool _graphShapesFull;

	float _zoomFactor;
	int _fluidShapeRadius;
	int _solidShapeRadius;

	Renderer(float zoomFactor = 1, float fluidSize = 1, float solidSize = 1, float searchRadius = 3);
	void init_solids(std::vector<Particle>* particles);
	void update_information(float time, float simTime, int numParticles, int numFluidParticles, float numUpdates,
		float avgDensity = -1, float cflNumber = -1, float avgSolverIters = -1., float watchedParticleDensity = -1, bool drawGraph = false);
	void update_arrows(Particle* watchedPartile);
	void draw(sf::RenderWindow* window, std::vector<Particle>* particles,
		int watchedParticleId, std::vector<int> markedParticlesId,
		std::vector<int> testedParticlesId, bool updateArrows = false, bool drawGraph = false, bool drawArrows = false);
};