// Paul Kull, 2022

#pragma once

#include <map>

#include "./Particle.h"

class HashManager {
  public:
	std::map<int, std::vector<Particle*>> _Buckets;
	int _hastableSize;
	int _prime1;
	int _prime2;
	int _cellSize;
	float _p1DIVcellSize;
	float _p2DIVcellSize;
	int _neighboringCells[9];

	HashManager(int radius = 10, int size = 3000, int prime1 = 73856093, int prime2 = 19349663);

	void insert_item(Particle* particle);
	int hash(sf::Vector2f pos);
	void reset_buckets();
	std::vector<Particle*> return_possible_neighbors(Particle* particle);
};