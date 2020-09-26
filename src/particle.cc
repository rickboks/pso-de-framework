#include "particle.h"
#include "particleupdatesettings.h"
#include "particleupdatemanager.h"
#include "util.h"
#include "rng.h"
#include <iostream>
#include <algorithm>
#include <functional>
#include <limits>
#include <random>
#include <cmath>
#include <IOHprofiler_experimenter.h>

Particle::Particle(int const D, ParticleUpdateSettings const*const settings)
	: Solution(D), v(D), p(D), g(D), pbest(std::numeric_limits<double>::max()), gbest(std::numeric_limits<double>::max()),
		settings(settings), psoCH(settings->psoCH){
	particleUpdateManager = updateManagers.at(settings->managerType)(x,v,p,g,settings->parameters,neighborhood);
}

Particle::Particle(Particle const & other)
	: Solution(other.D), v(other.v), p(other.p), g(other.g),
	pbest(other.pbest), gbest(other.gbest), 
	neighborhood(other.neighborhood), particleUpdateManager(NULL),
	settings(other.settings), psoCH(other.psoCH){

	x = other.x;
	evaluated=other.evaluated;
	fitness = other.fitness;

	if (settings != NULL) // This constructor is used also by DE
		particleUpdateManager = updateManagers.at(settings->managerType)(x,v,p,g,settings->parameters,neighborhood);
}

Particle::~Particle(){
	if (particleUpdateManager != NULL)
		delete particleUpdateManager;
}

std::vector<double> Particle::getV() const {
	return v;
}

double Particle::getV(int const dim) const {
	return v[dim];
}

void Particle::setV(std::vector<double> v){
	this->v = v;
}

void Particle::setV(int const dim, double val){
	this->v[dim] = val;
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
	std::vector<double> const oldV = v;
	std::vector<double> const oldX = x;
	int resamples = 0;

	while(true){
		particleUpdateManager->updateVelocity(progress);
		psoCH->repairVelocityPre(this);
		particleUpdateManager->updatePosition();
		if (psoCH->resample(this, resamples)){
			x = oldX; // reset position and velocity
			v = oldV;
			resamples++;
		} else 
			break;
	}
	psoCH->repair(this); // Generic repair
}

double Particle::getGbest() const {
	return gbest;
}

double Particle::getPbest() const {
	return pbest;
}

double Particle::getP(int const i) const {
	return p[i];
}

std::vector<double> Particle::getG() const {
	return g;
}

std::vector<double> Particle::getP() const {
	return p;
}

void Particle::updateGbest(){
	int bestNeighbor = -1;

	if (fitness < gbest){ // First check own fitness
		gbest = fitness;
		g = x;
	}

	double bestScore = gbest;
	for (unsigned int i = 0; i < neighborhood.size(); i++){ // Check neighbors fitness
		double const currentScore = neighborhood[i]->getPbest();
		if (currentScore < bestScore){
			bestScore = currentScore;
			bestNeighbor = i;
		}
	}

	if (bestNeighbor != -1){ // Update gbest
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

int Particle::getNumberOfNeighbors() const{
	return neighborhood.size();
}

void Particle::setXandUpdateV(std::vector<double> x, double fitness){
	subtract(x, this->x, v); // Reverse engineer velocity
	this->x = x;
	this->fitness = fitness;
}
