// Paul Kull, 2022

#pragma once

#include <map>

#include "./Functions.h"
#include "./Particle.h"

class CompactHashManager {
public:
	std::map<std::pair<unsigned int, unsigned int>, std::vector<Particle*>*> _map;
	std::vector<std::vector<Particle*>> _buckets;
	int _numParticles;
	int _areaWidth;
	int _areaHeight;
	int _prime1;
	int _prime2;
	int _cellSize;
	float _p1DIVcellSize;
	float _p2DIVcellSize;
	float _distanceNorm;
	int _index;
	int _bucketCounter;
	std::pair<unsigned int, unsigned int> _particleCell;
	sf::Vector2f _distance;
	std::pair<unsigned int, unsigned int> _neighboringCells[9];

	CompactHashManager(int radius = 10, int areaWidth = 1, int areaHeight = 1);


	void insert_item(Particle* particle);
	std::pair<unsigned int, unsigned int> hash(sf::Vector2f pos);
	void reset_buckets();
	void reset_map();
	std::vector<Particle*> return_neighbors(Particle* particle, float radius);
	void update(int numParticles);
};