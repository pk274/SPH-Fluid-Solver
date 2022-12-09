// Paul Kull, 2022

#pragma once

#include <map>

#include "./Functions.h"
#include "./Particle.h"

class HashManager {
  public:
	std::map<std::pair<unsigned int, unsigned int>, std::vector<Particle*>> _Buckets;
	int _areaWidth;
	int _areaHeight;
	int _prime1;
	int _prime2;
	int _cellSize;
	float _p1DIVcellSize;
	float _p2DIVcellSize;
	float _distanceNorm;
	int _index;
	std::pair<unsigned int, unsigned int> _particleCell;
	sf::Vector2f _distance;
	std::pair<unsigned int, unsigned int> _neighboringCells[9];

	HashManager(int radius = 10, int areaWidth = 1, int areaHeight = 1);

	void insert_item(Particle* particle);
	std::pair<unsigned int, unsigned int> hash(sf::Vector2f pos);
	void reset_buckets();
	std::vector<Particle*> return_neighbors(Particle* particle, float radius);
};