// Paul Kull, 2022

#include <iostream>
#include "CompactHashManager.h"

// __________________________________________________________________________
CompactHashManager::CompactHashManager(int radius, int areaWidth, int areaHeight) {
	_cellSize = radius;
	_areaWidth = areaWidth;
	_areaHeight = areaHeight;
	_map = std::map<std::pair<unsigned int, unsigned int>, std::vector<Particle*>*>();
	for (int i = 0; i < _areaWidth; i++) {
		for (int j = 0; j < _areaHeight; j++) {
			_map.insert(std::pair<std::pair<unsigned int, unsigned int>,
				std::vector<Particle*>*>(std::pair(i, j), NULL));
		}
	}
	_buckets = std::vector<std::vector<Particle*>>();
	reset_map();
	reset_buckets();
}


// __________________________________________________________________________
void CompactHashManager::insert_item(Particle* particle) {
	std::pair<unsigned int, unsigned int> bucketHash = hash(particle->_position);
	if (_map.at(bucketHash) == NULL) {
		_map.at(bucketHash) = &_buckets[_bucketCounter];
		_map.at(bucketHash)->push_back(particle);
		_bucketCounter += 1;
	}
	else {
		_map.at(bucketHash)->push_back(particle);
	}
	if (_bucketCounter > _numParticles) {
		throw std::runtime_error("bucketCounter > numParticles, please contact Simon G. for help");
	}
}

// __________________________________________________________________________
void CompactHashManager::reset_map() {
	for (int i = 0; i < _areaWidth; i++) {
		for (int j = 0; j < _areaHeight; j++) {
			_map.at(std::pair(i, j)) = NULL;
		}
	}
}

// __________________________________________________________________________
void CompactHashManager::reset_buckets() {
	_buckets.clear();
	for (int i = 0; i < _numParticles; i++) {
		_buckets.push_back(std::vector<Particle*>());
	}
}

// __________________________________________________________________________
std::pair<unsigned int, unsigned int> CompactHashManager::hash(sf::Vector2f pos) {

	std::pair<unsigned int, unsigned int> hashvalue = std::pair(((unsigned int)(std::floor(pos.x / _cellSize)) % _areaWidth),
		((unsigned int)(std::floor(pos.y / _cellSize)) % _areaHeight));

	return hashvalue;
}


// __________________________________________________________________________
std::vector<Particle*> CompactHashManager::return_neighbors(Particle* particle, float radius) {
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
		if (!_map.contains(_neighboringCells[i])) { continue; }
		Bucket = _map.find(_neighboringCells[i])->second;
		if (Bucket == NULL) { continue; }
		for (int j = 0; j < Bucket->size(); j++) {
			_distance = Functions::calculate_distance(Bucket->at(j)->_position, particle->_position);
			_distanceNorm = Functions::calculate_distance_norm(_distance);
			if (_distanceNorm < radius) {
				neighbors.push_back(Bucket->at(j));
			}
		}
	}

	return neighbors;
}

// __________________________________________________________________________
void CompactHashManager::update(int numParticles) {
	if (numParticles != _numParticles) {
		_numParticles = numParticles;
	}
	_bucketCounter = 0;
	reset_buckets();
	reset_map();

}