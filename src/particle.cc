#include "particle.h"
#include "particleupdatesettings.h"
#include "particleupdatemanager.h"
#include "rng.h"
#include <iostream>
#include <algorithm>
#include <functional>
#include <random>
#include <cmath>

constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();
constexpr double DOUBLE_MIN = std::numeric_limits<double>::min();

Particle::Particle(int const D, ParticleUpdateSettings& particleUpdateSettings)
	: D(D) , x(D), v(D), p(D), pbest(DOUBLE_MAX), g(D), gbest(DOUBLE_MAX), settings(particleUpdateSettings),
	vMax(particleUpdateSettings.vMax), evaluated(false){

	particleUpdateManager = ParticleUpdateManagerFactory::createParticleUpdateManager(x,v,p,g,particleUpdateSettings,neighborhood);
}

Particle::Particle(const Particle& other)
: D(other.D), x(other.x), v(other.v), p(other.p),
fitness(other.fitness), pbest(other.pbest), g(other.g), gbest(other.gbest), neighborhood(other.neighborhood),
settings(other.settings), vMax(other.vMax), evaluated(other.evaluated){
	particleUpdateManager = ParticleUpdateManagerFactory::createParticleUpdateManager(x,v,p,g,settings,neighborhood);
}


Particle::~Particle(){
	delete particleUpdateManager;
}

void Particle::randomize(std::vector<double> lowerBounds, std::vector<double> upperBounds){
	for (int i = 0; i < D; i++){
		x[i] = rng.randDouble(lowerBounds[i], upperBounds[i]);
		v[i] = 0;
	}
	p = x;
	g = x;

	evaluated = false;
}

std::vector<double> Particle::getPosition() const {
	return x;
}

void Particle::setPosition(std::vector<double> position, double fitness){
	evaluated = true;
	this->fitness = fitness;
	this->x = position;
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

double Particle::getResultingVelocity() const{
	double velocity = 0;

	for (int i = 0; i < (int)v.size(); i++){
		velocity += v[i] * v[i];
	}
	velocity = sqrt(velocity);
	return velocity;
}

double Particle::getFitness() const{
	return fitness;
}

double Particle::evaluate(evaluate_function_t evalFunc){
	if (!evaluated){
		evaluated = true;
		double f;
		evalFunc(&x[0], &f);
		this->fitness = f;
		return f;
	} else {
		return fitness;
	}
}

void Particle::replaceNeighbors(std::map<Particle*, Particle*> mapping){
	for (int i = 0; i < (int) neighborhood.size(); i++){
		neighborhood[i] = mapping[neighborhood[i]];
	}
}