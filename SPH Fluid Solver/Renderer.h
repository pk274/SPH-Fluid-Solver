#pragma once
// Paul Kull, 2022

#include <SFML/Graphics.hpp>
#include "./SolidParticle.h"
#include "./FluidParticle.h"

class Renderer {
  public:
	std::vector<sf::CircleShape> _particleShapes;
	sf::CircleShape _searchRadiusShape;
	sf::Vector2f _searchRadiusOffset;
	sf::Font _font;
	sf::Text _description;
	sf::Text _information;

	float _zoomFactor;
	int _fluidShapeRadius;
	int _solidShapeRadius;

	Renderer(float zoomFactor = 1, int fluidRadius = 1, int solidRadius = 1, int searchRadius = 3);
	void update_graphics(std::vector<Particle>* particles, int watchedParticleId, std::vector<int> markedParticlesId);
	void update_information(int numParticles, int numUpdates);
	void draw(sf::RenderWindow* window);
};