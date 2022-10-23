// Paul Kull, 2022

#include "HashManager.h"

// __________________________________________________________________________
HashManager::HashManager(int radius, int hashtableSize, int prime1, int prime2) {
	_cellSize = radius;
	_hastableSize = hashtableSize;
	_prime1 = prime1;
	_prime2 = prime2;
	_p1DIVcellSize = _prime1 / _cellSize;
	_p2DIVcellSize = _prime2 / _cellSize;
	_Buckets = std::map<int, std::vector<Particle*>>();
	reset_buckets();
}


// __________________________________________________________________________
void HashManager::insert_item(Particle* particle) {
	int bucketHash = hash(particle->_position);
	_Buckets[bucketHash].push_back(particle);
}


// __________________________________________________________________________
void HashManager::reset_buckets() {
	_Buckets.clear();
	for (int i = 0; i < _hastableSize; i++) {
		_Buckets.insert(std::pair(i, std::vector<Particle*>()));
	}
}

// __________________________________________________________________________
int HashManager::hash(sf::Vector2f pos) {
	if (pos.x < 0) { pos.x = 0; }
	if (pos.y < 0) { pos.y = 0; }
	return ((int)(std::floor(pos.x / _cellSize) * _prime1) ^
		(int)(std::floor(pos.y / _cellSize) * _prime2)) % _hastableSize;
}

// __________________________________________________________________________
std::vector<Particle*> HashManager::return_possible_neighbors(Particle* particle) {
	int index = 0;
	for (int i = -_cellSize; i <= _cellSize; i += _cellSize) {
		for (int j = -_cellSize; j <= _cellSize; j += _cellSize) {
			int hashvalue = hash(particle->_position + sf::Vector2f(i, j));
			_neighboringCells[index] = hashvalue;
			index += 1;
		}
	}
	std::vector<Particle*> possibleNeighbors = std::vector<Particle*>();
	for (int i = 0; i < 9; i++) {
		std::vector<Particle*>* Bucket = &_Buckets.find(_neighboringCells[i])->second;
		for (int j = 0; j < Bucket->size(); j++) {
			possibleNeighbors.push_back(Bucket->at(j));
		}
	}
	return possibleNeighbors;
}