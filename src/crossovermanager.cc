#include "crossovermanager.h"
#include "genome.h"
#include "rng.h"
#include <iostream>
#include <algorithm>

// FACTORY
CrossoverManager* CrossoverManagerFactory::createCrossoverManager(CrossoverType const crossoverType, 
	std::vector<Genome*>& genomes, std::vector<Genome*>& mutants,double const Cr, int const D){
	switch(crossoverType){
		case BINOMIAL:
			return new BinomialCrossoverManager(genomes, mutants, Cr, D);
		case EXPONENTIAL:
			return new ExponentialCrossoverManager(genomes, mutants, Cr, D);
		default:
			throw std::invalid_argument("Error: Invalid DE crossover type");
	}
}

// BASE
CrossoverManager::CrossoverManager(std::vector<Genome*>& genomes, std::vector<Genome*>& mutants, 
	double const Cr, int const D)
: genomes(genomes), mutants(mutants), Cr(Cr), D(D){

}

CrossoverManager::~CrossoverManager(){

}

// BINOMIAL

BinomialCrossoverManager::BinomialCrossoverManager(std::vector<Genome*>& genomes, std::vector<Genome*>& mutants, 
	double const Cr, int const D)
: CrossoverManager(genomes, mutants, Cr, D){

}

std::vector<Genome*> BinomialCrossoverManager::crossover() {
	std::vector<Genome*> donors;
	donors.reserve(genomes.size());
	std::vector<double> x(D);

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<double> mutantX = mutants[i]->getPosition();
		std::vector<double> parentX = genomes[i]->getPosition();
		int const jrand = rng.randInt(0,D-1);
		for (int j = 0; j < D; j++){
			if (j == jrand || rng.randDouble(0,1) < Cr){
				x[j] = mutantX[j];
			} else {
				x[j] = parentX[j]; 
			}
		}

		Genome* donor = new Genome(D);
		donor->setPosition(x);
		donors.push_back(donor);
	}

	return donors;
}

// Exponential
ExponentialCrossoverManager::ExponentialCrossoverManager(std::vector<Genome*>& genomes, std::vector<Genome*>& mutants, 
	double const Cr, int const D)
: CrossoverManager(genomes, mutants, Cr, D){

}

std::vector<Genome*> ExponentialCrossoverManager::crossover() {
	std::vector<Genome*> donors;
	donors.reserve(genomes.size());
	std::vector<double> x(D);

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<double> mutantX = mutants[i]->getPosition();
		std::vector<double> parentX = genomes[i]->getPosition();
		std::vector<int> mutantIndices;
		int const n = rng.randInt(0,D-1);

		int L = 0;
		do {
			L++;
		} while (rng.randDouble(0,1) < Cr && L < D);

		int steps = 0;
		for (int i = n; steps < L; i++){
			mutantIndices.push_back(i%D);
			steps++;
		}

		for (int i = 0; i < D; i++){
			if (std::find(mutantIndices.begin(), mutantIndices.end(), i) != mutantIndices.end())
				x[i] = mutantX[i];
			else 
				x[i] = parentX[i];
		}

		Genome* donor = new Genome(D);
		donor->setPosition(x);
		donors.push_back(donor);
	}

	return donors;
}