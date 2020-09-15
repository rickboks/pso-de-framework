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

Particle::Particle(int const D, ParticleUpdateSettings& particleUpdateSettings)
	: x(D), v(D), p(D), g(D), pbest(std::numeric_limits<double>::max()), gbest(std::numeric_limits<double>::max()), evaluated(false),
		settings(particleUpdateSettings), psoCH(settings.psoCH), D(D), isPSO(true){

	particleUpdateManager = ParticleUpdateManager::createParticleUpdateManager(x,v,p,g,particleUpdateSettings,neighborhood);
}

Particle::Particle(const Particle& other)
	: x(other.x), v(other.v), p(other.p), g(other.g),
	pbest(other.pbest), gbest(other.gbest), evaluated(other.evaluated), 
	fitness(other.fitness), neighborhood(other.neighborhood),
	settings(other.settings), psoCH(other.psoCH), D(other.D), isPSO(other.isPSO){

	particleUpdateManager = ParticleUpdateManager::createParticleUpdateManager(x,v,p,g,settings,neighborhood);
}

//for DE
Particle::Particle(int const D)
	: x(D), evaluated(false), fitness(std::numeric_limits<double>::max()), particleUpdateManager(NULL), psoCH(NULL), D(D),  isPSO(false){
}

//When using this construtor, note that only the position is initialized
Particle::Particle(std::vector<double> x)
	: x(x), evaluated(false), fitness(std::numeric_limits<double>::max()), particleUpdateManager(NULL),psoCH(NULL),D(x.size()),isPSO(false){
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

void Particle::setX(std::vector<double> x, double fitness, bool updateVelocity){
	if (updateVelocity)
		subtract(x, this->x, v); // Reverse engineer velocity
	this->x = x;
	this->fitness = fitness;
}

void Particle::setX(std::vector<double> x){
	this->x = x;
}

double Particle::getFitness() const{
	return fitness;
}

void Particle::setFitness(double const f){
	this->fitness = f;
}

double Particle::evaluate(std::shared_ptr<IOHprofiler_problem<double> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
	if (!evaluated){
		evaluated = true;		
		fitness = problem->evaluate(x);
		logger->do_log(problem->loggerCOCOInfo());
		return fitness;
	} else {
		return fitness;
	}
}

std::vector<double> Particle::getX() const {
	return x;
}

std::string Particle::positionString() const {
	std::string pos = "";
	for (int i = 0; i < D -1; i++){
		pos += std::to_string(x[i]);
		pos += " ";
	}
	pos += std::to_string(x[D-1]);

	return pos;
}

void Particle::randomize(std::vector<double> const lowerBounds, std::vector<double> const upperBounds){
	for (int i = 0; i < D; i++){
		x[i] = rng.randDouble(lowerBounds[i], upperBounds[i]);
	}
	evaluated=false;
}

bool Particle::operator < (const Particle& s) const {
	return fitness < s.getFitness();
}

void Particle::setX(int const dim, double const val){
	x[dim] = val;
}

double Particle::getX(int const dim) const {
	return x[dim];
}
