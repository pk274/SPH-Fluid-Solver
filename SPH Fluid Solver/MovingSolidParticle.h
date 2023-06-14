// Paul Kull, 2023
#pragma once


#include <SFML/Graphics.hpp>

#include "./Particle.h"

class MovingSolidParticle : public Particle {
public:
	static const float _size;
	static const float _mass;
	std::vector<sf::Vector2f> _velocities;
	std::vector<sf::Vector2f> _triggers;
	std::vector<Particle> _friends;
	int _statusPointer;

	MovingSolidParticle(int id, sf::Vector2f pos, sf::Vector2f vel,
		std::vector<sf::Vector2f> velocities, std::vector<sf::Vector2f> triggers, std::vector<int> friends);
};