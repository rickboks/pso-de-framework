#include "crossovermanager.h"
#include<unordered_set>

CrossoverManager::CrossoverManager(int const D): D(D){}
CrossoverManager::~CrossoverManager(){}

#define LC(X) [](int const D){return new X(D);}
std::map<std::string, std::function<CrossoverManager* (int const)>> const crossovers({
		{"B", LC(BinomialCrossoverManager)},
		{"E", LC(ExponentialCrossoverManager)},
});

BinomialCrossoverManager::BinomialCrossoverManager(int const D) : CrossoverManager(D){}

std::vector<Particle*> BinomialCrossoverManager::crossover(std::vector<Particle*>const& genomes, std::vector<Particle*>const& mutants, std::vector<double>const& Crs){
	std::vector<Particle*> trials;
	trials.reserve(genomes.size());
	std::vector<double> x(D);

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<double> const mutantX = mutants[i]->getX();
		std::vector<double> const parentX = genomes[i]->getX();
		trials.push_back(new Particle(singleCrossover(genomes[i]->getX(), mutants[i]->getX(), Crs[i])));
	}

	return trials;
}

std::vector<double> BinomialCrossoverManager::singleCrossover(std::vector<double>const& target, std::vector<double>const& donor, double const Cr){
	std::vector<double> x(D);
	int const jrand = rng.randInt(0,D-1);
	for (int j = 0; j < D; j++){
		if (j == jrand || rng.randDouble(0,1) < Cr){
			x[j] = donor[j];
		} else {
			x[j] = target[j]; 
		}
	}
	return x;
}

ExponentialCrossoverManager::ExponentialCrossoverManager(int const D): CrossoverManager(D){}

std::vector<Particle*> ExponentialCrossoverManager::crossover(std::vector<Particle*>const& genomes, std::vector<Particle*>const& mutants, std::vector<double>const& Crs){
	std::vector<Particle*> trials;
	for (unsigned int i = 0; i < genomes.size(); i++){
		trials.push_back(new Particle(singleCrossover(genomes[i]->getX(), mutants[i]->getX(), Crs[i])));
	}
	return trials;
}

std::vector<double> ExponentialCrossoverManager::singleCrossover(std::vector<double>const& target, std::vector<double>const& donor, double const Cr){
	std::vector<double> x(D);
	int const start = rng.randInt(0,D-1);

	int L = 0;
	do {
		L++;
	} while (rng.randDouble(0,1) <= Cr && L <= D);

	int const end = (start+L-1) % D;

	for (int i = 0; i < D; i++){
		if ((end >= start && (i >= start && i <= end)) || (start > end && (i <= end || i >= start)))
			x[i] = donor[i];
		else 
			x[i] = target[i];
	}

	return x;
}
