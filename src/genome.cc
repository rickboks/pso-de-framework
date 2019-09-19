#include "genome.h"
#include "rng.h"
#include <random>
#include <algorithm>
#include <functional>
#include "particle.h"
#include <iostream>

Genome::Genome(int const D): 
	D(D), x(D), evaluated(false){

}

Genome::Genome(Particle* particle) :
	D(particle->getPosition().size()), x(particle->getPosition()), evaluated(true), fitness(particle->getFitness()){

}

Genome::Genome(std::vector<double>x): 
	D(x.size()), x(x), evaluated(false){

}

void Genome::randomize(std::vector<double> lowerBounds, std::vector<double> upperBounds){
	for (int i = 0; i < D; i++){
		x[i] = rng.randDouble(lowerBounds[i], upperBounds[i]);
	}
}

double Genome::getFitness() const {
	return fitness;
}

std::vector<double> Genome::getPosition() const {
	return x;
}

void Genome::setPosition(std::vector<double> const x){
	this->x = x;
}

int Genome::getDimension() const {
	return D;
}

void Genome::print() const {
	for (double d : x){
		std::cout << d << " ";
	}

	std::cout << std::endl;
}


double Genome::evaluate(evaluate_function_t evalFunc){
	if (!evaluated){
		evaluated = true;
		double fitness;
		evalFunc(&x[0], &fitness);
		this->fitness = fitness;
		return fitness;
	} else {
		return fitness;
	}
}

bool Genome::operator < (const Genome& genome) const
{
	return fitness < genome.getFitness();
}

void Genome::setWeightFactor(double weightFactor){
	this->weightFactor = weightFactor;
}

double Genome::getWeightFactor() const{
	return weightFactor;
}