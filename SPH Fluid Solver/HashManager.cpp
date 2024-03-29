// Paul Kull, 2022

#include <iostream>
#include "HashManager.h"

// __________________________________________________________________________
HashManager::HashManager(int radius, int areaWidth, int areaHeight) {
	_cellSize = radius;
	_areaWidth = areaWidth;
	_areaHeight = areaHeight;
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
	for (int i = 0; i < _areaWidth; i++) {
		for (int j = 0; j < _areaHeight; j++) {
			_Buckets.insert(std::pair(std::pair<int, int>(i, j), std::vector<Particle*>()));
		}
	}
}

// __________________________________________________________________________
std::pair<unsigned int, unsigned int> HashManager::hash(sf::Vector2f pos) {

	std::pair<unsigned int, unsigned int> hashvalue = std::pair(((unsigned int)(std::floor(pos.x / _cellSize)) % _areaWidth),
		((unsigned int)(std::floor(pos.y / _cellSize)) % _areaHeight));

	return hashvalue;
}


// __________________________________________________________________________
std::vector<Particle*> HashManager::return_neighbors(Particle* particle, float radius) {
	_particleCell = hash(particle->_position);
	if (!(_particleCell.first == _neighboringCells[4].first && _particleCell.second == _neighboringCells[4].second)) {
		_index = 0;
		for (int i = -_cellSize; i <= _cellSize; i += _cellSize) {
			for (int j = -_cellSize; j <= _cellSize; j += _cellSize) {
				std::pair<unsigned int, unsigned int> hashvalue = hash(particle->_position + sf::Vector2f(i, j));
				_neighboringCells[_index] = hashvalue;
				_index += 1;
			}
		}
	}

	std::vector<Particle*> neighbors = std::vector<Particle*>();
	std::vector<Particle*>* Bucket;

	for (int i = 0; i < 9; i++) {
		if (!_Buckets.contains(_neighboringCells[i])) { continue; }
		Bucket = &_Buckets.find(_neighboringCells[i])->second;
		for (int j = 0; j < Bucket->size(); j++) {
			_distance = Functions::calculate_distance(Bucket->at(j)->_position, particle->_position);
			_distanceNorm = Functions::calculate_distance_norm(_distance);
			if (_distanceNorm < radius) {
				neighbors.push_back(Bucket->at(j));
			}
		}
	}
	std::sort(neighbors.begin(), neighbors.end());
	neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());

	return neighbors;
}
