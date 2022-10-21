#pragma once
// Paul Kull, 2022

#include <SFML/Graphics.hpp>
#include "./SolidParticle.h"
#include "./FluidParticle.h"

class Renderer {
  public:
	std::vector<sf::CircleShape> _particleShapes;
	int _fluidParticleRadius;
	int _solidParticleRadius;

	Renderer(int fluidRadius = 5, int solidRadius = 5);
	void update_graphics(std::vector<Particle>* particles);
	void draw(sf::RenderWindow* window);
};