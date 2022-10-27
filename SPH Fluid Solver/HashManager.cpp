// Paul Kull, 2022

#include <iostream>
#include "HashManager.h"

// __________________________________________________________________________
HashManager::HashManager(int radius, int hashtableSize, int prime1, int prime2) {
	_cellSize = radius;
	_hastableSize = hashtableSize;
	_sqrtHashtableSize = std::ceil(std::sqrt(hashtableSize));
	_prime1 = prime1;
	_prime2 = prime2;
	_p1DIVcellSize = _prime1 / _cellSize;
	_p2DIVcellSize = _prime2 / _cellSize;
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

// __________________________________________________________________________________
sf::Vector2f calculate_distance(sf::Vector2f pos1, sf::Vector2f pos2) {
	return sf::Vector2f(pos1.x - pos2.x, pos1.y - pos2.y);
}


// __________________________________________________________________________________
float calculate_distance_norm(sf::Vector2f distance) {
	return std::abs(std::sqrt(std::pow(distance.x, 2) + std::pow(distance.y, 2)));
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
			distance = calculate_distance(Bucket->at(j)->_position, particle->_position);
			distanceNorm = calculate_distance_norm(distance);
			if (distanceNorm <= radius) {
				Bucket->at(j)->_distance = distance;
				Bucket->at(j)->_distanceNorm = distanceNorm;
				neighbors.push_back(Bucket->at(j));
			}
		}
	}
	return neighbors;
}