// Paul Kull, 2022

#include <iostream>
#include "HashManager.h"

// __________________________________________________________________________
HashManager::HashManager(int radius, int hashtableSize) {
	_cellSize = radius * 2;
	_hastableSize = hashtableSize;
	_sqrtHashtableSize = std::ceil(std::sqrt(hashtableSize));
	_Buckets = std::map<std::pair<unsigned int, unsigned int>, std::vector<Particle*>>();
	reset_buckets();
}


// __________________________________________________________________________
void HashManager::insert_item(Particle* particle) {
	std::pair<unsigned int, unsigned int> bucketHash = hash(particle->_position);
	_Buckets[bucketHash].push_back(particle);
}


// __________________________________________________________________________
void HashManager::reset_buckets() {
	_Buckets.clear();
	for (int i = 0; i < _sqrtHashtableSize; i++) {
		for (int j = 0; j < _sqrtHashtableSize; j++) {
			_Buckets.insert(std::pair(std::pair<int, int>(i, j), std::vector<Particle*>()));
		}
	}
}

// __________________________________________________________________________
std::pair<unsigned int, unsigned int> HashManager::hash(sf::Vector2f pos) {

	std::pair<unsigned int, unsigned int> hashvalue = std::pair(((unsigned int)(std::floor(pos.x / _cellSize)) % _sqrtHashtableSize),
		((unsigned int)(std::floor(pos.y / _cellSize)) % _sqrtHashtableSize));

	return hashvalue;
}


// __________________________________________________________________________
std::vector<Particle*> HashManager::return_neighbors(Particle* particle, float radius) {
	int index = 0;
	for (int i = -_cellSize; i <= _cellSize; i += _cellSize) {
		for (int j = -_cellSize; j <= _cellSize; j += _cellSize) {
			std::pair<unsigned int, unsigned int> hashvalue = hash(particle->_position + sf::Vector2f(i, j));
			_neighboringCells[index] = hashvalue;
			index += 1;
		}
	}

	std::vector<Particle*> neighbors = std::vector<Particle*>();
	std::vector<Particle*>* Bucket;
	sf::Vector2f distance;
	float distanceNorm;
	for (int i = 0; i < 9; i++) {
		if (!_Buckets.contains(_neighboringCells[i])) { continue; }
		Bucket = &_Buckets.find(_neighboringCells[i])->second;
		for (int j = 0; j < Bucket->size(); j++) {
			if (Bucket->at(j)->_id == particle->_id) { continue; }
			distance = Functions::calculate_distance(Bucket->at(j)->_position, particle->_position);
			distanceNorm = Functions::calculate_distance_norm(distance);
			if (distanceNorm <= radius) {
				neighbors.push_back(Bucket->at(j));
			}
		}
	}
	std::sort(neighbors.begin(), neighbors.end());
	neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());

	return neighbors;
}