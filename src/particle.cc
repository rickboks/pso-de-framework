#include "particle.h"
#include "particleupdatesettings.h"
#include "particleupdatemanager.h"
#include "rng.h"
#include <iostream>
#include <algorithm>
#include <functional>
#include <random>
#include <cmath>


Particle::Particle(int const D, ParticleUpdateSettings& particleUpdateSettings)
	: Solution(D), v(D), p(D), pbest(std::numeric_limits<double>::max()), g(D), gbest(std::numeric_limits<double>::max()), 
		settings(particleUpdateSettings), vMax(particleUpdateSettings.vMax){

	particleUpdateManager = ParticleUpdateManagerFactory::createParticleUpdateManager(x,v,p,g,particleUpdateSettings,neighborhood);
}

Particle::Particle(const Particle& other)
: 	Solution(other.x, other.fitness), v(other.v), p(other.p),
	pbest(other.pbest), g(other.g), gbest(other.gbest), neighborhood(other.neighborhood),
	settings(other.settings), vMax(other.vMax){
	particleUpdateManager = ParticleUpdateManagerFactory::createParticleUpdateManager(x,v,p,g,settings,neighborhood);
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
	double bestScore = gbest;
	int const size = neighborhood.size();

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

		if (fitness < gbest){
			gbest = fitness;
			g = x;
		}
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
