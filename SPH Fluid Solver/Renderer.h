#pragma once
// Paul Kull, 2022

#include <SFML/Graphics.hpp>
#include "./SolidParticle.h"
#include "./FluidParticle.h"
#include "./Functions.h"

class Renderer {
  public:
	std::vector<sf::CircleShape> _particleShapes;
	std::vector<sf::RectangleShape>_arrowBodies;
	std::vector<sf::CircleShape> _arrowHeads;
	sf::CircleShape _searchRadiusShape;
	sf::Vector2f _searchRadiusOffset;
	sf::Font _font;
	sf::Text _description;
	sf::Text _information;

	float _zoomFactor;
	int _fluidShapeRadius;
	int _solidShapeRadius;

	Renderer(float zoomFactor = 1, float fluidSize = 1, float solidSize = 1, float searchRadius = 3);
	void update_graphics(std::vector<Particle>* particles, int watchedParticleId, std::vector<int> markedParticlesId, std::vector<int> testedParticlesId);
	void update_information(int numParticles, float numUpdates);
	void update_arrows(std::vector<Particle>* particles, Particle* watchedPartile);
	void draw(sf::RenderWindow* window);
};