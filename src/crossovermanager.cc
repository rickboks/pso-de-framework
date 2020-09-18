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
	std::vector<double> x(this->D);

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<double> mutantX = mutants[i]->getX();
		std::vector<double> parentX = genomes[i]->getX();
		trials.push_back(new Particle(singleCrossover(genomes[i]->getX(), mutants[i]->getX(), Crs[i])));
	}

	return trials;
}

std::vector<double> BinomialCrossoverManager::singleCrossover(std::vector<double>const& target, std::vector<double>const& donor, double const Cr){
	std::vector<double> x(this->D);
	int const jrand = rng.randInt(0,this->D-1);
	for (int j = 0; j < this->D; j++){
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
	std::vector<double> x(this->D);
	std::unordered_set<int> mutantIndices;
	int const n = rng.randInt(0,this->D-1);

	int L = 0;
	do {
		L++;
	} while (rng.randDouble(0,1) < Cr && L < this->D);

	int steps = 0;
	for (int i = n; steps < L; i++){
		mutantIndices.insert(i%this->D);
		steps++;
	}

	for (int i = 0; i < this->D; i++) {
		if (mutantIndices.find(i) != mutantIndices.end()){
			x[i] = donor[i];
		} else {
			x[i] = target[i];
		}
	}

	return x;
}
