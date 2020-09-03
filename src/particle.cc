#include "particle.h"
#include "particleupdatesettings.h"
#include "particleupdatemanager.h"
#include "util.h"
#include "rng.h"
#include <iostream>
#include <algorithm>
#include <functional>
#include <random>
#include <cmath>

Particle::Particle(int const D, ParticleUpdateSettings& particleUpdateSettings)
	: Solution(D), v(D), p(D), g(D), pbest(std::numeric_limits<double>::max()), gbest(std::numeric_limits<double>::max()), 
		settings(particleUpdateSettings), vMax(particleUpdateSettings.vMax){
	particleUpdateManager = ParticleUpdateManager::createParticleUpdateManager(x,v,p,g,particleUpdateSettings,neighborhood);
}

Particle::Particle(const Particle& other)
: 	Solution(other.x, other.fitness), v(other.v), p(other.p), g(other.g),
	pbest(other.pbest),  gbest(other.gbest), neighborhood(other.neighborhood),
	settings(other.settings), vMax(other.vMax){
	particleUpdateManager = ParticleUpdateManager::createParticleUpdateManager(x,v,p,g,settings,neighborhood);
}

//When using this construtor, note that only the position is initialized
Particle::Particle(std::vector<double> x)
: 	Solution(x), particleUpdateManager(NULL) {
}

Particle::~Particle(){
	if (particleUpdateManager != NULL)
		delete particleUpdateManager;
}

std::vector<double> Particle::getVelocity() const {
	return v;
}

void Particle::setVelocity(std::vector<double> v){
	this->v = v;
}

void Particle::addNeighbor(Particle* const neighbor){
	neighborhood.push_back(neighbor);
}

void Particle::removeNeighbor(Particle* const neighbor){
	neighborhood.erase(std::remove(neighborhood.begin(), neighborhood.end(), neighbor), neighborhood.end());
}

void Particle::removeAllNeighbors(){
	neighborhood.clear();
}

void Particle::updateVelocityAndPosition(double progress){
	evaluated = false;
	particleUpdateManager->update(progress);
}

double Particle::getGbest() const {
	return gbest;
}

double Particle::getPbest() const {
	return pbest;
}

std::vector<double> Particle::getG() const {
	return g;
}

std::vector<double> Particle::getP() const {
	return p;
}

void Particle::updateGbest(){
	int bestNeighbor = -1;
	int const size = neighborhood.size();

	if (fitness < gbest){
		gbest = fitness;
		g = x;
	}

	double bestScore = gbest;
	for (int i = 0; i < size; i++){
		double const currentScore = neighborhood[i]->getPbest();
		if (currentScore < bestScore){
			bestScore = currentScore;
			bestNeighbor = i;
		}
	}

	if (bestNeighbor != -1){
		gbest = bestScore;
		g = neighborhood[bestNeighbor]->getG();
	}
}

void Particle::updatePbest(){
	if (fitness < pbest){
		pbest = fitness;
		p = x;
	}
}

bool Particle::isNeighbor(Particle* particle) const {
	return std::find(neighborhood.begin(), neighborhood.end(), particle) != neighborhood.end();
}

void Particle::replaceNeighbors(std::map<Particle*, Particle*> mapping){
	for (int i = 0; i < (int) neighborhood.size(); i++){
		neighborhood[i] = mapping[neighborhood[i]];
	}
}

int Particle::getAmountOfNeighbors(){
	return neighborhood.size();
}

void Particle::setPosition(std::vector<double> x, double fitness, bool updateVelocity){
	if (updateVelocity)
		subtract(x, this->x, v); // Reverse engineer velocity
	this->x = x;
	this->fitness = fitness;
}
