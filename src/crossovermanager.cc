#include "util.h"
#include "crossovermanager.h"
#include<unordered_set>

CrossoverManager::CrossoverManager(int const D): D(D){}
CrossoverManager::~CrossoverManager(){}

std::vector<Solution*> CrossoverManager::crossover(std::vector<Solution*>const& genomes, std::vector<Solution*>const& mutants, std::vector<double>const& Crs) const{
	std::vector<Solution*> trials;
	trials.reserve(genomes.size());
	std::vector<double> x(D);

	for (unsigned int i = 0; i < genomes.size(); i++){
		trials.push_back(new Solution(singleCrossover(genomes[i]->getX(), mutants[i]->getX(), Crs[i])));
	}
	return trials;
}

#define LC(X) [](int const D){return new X(D);}
std::map<std::string, std::function<CrossoverManager* (int const)>> const crossovers({
		{"B", LC(BinomialCrossoverManager)},
		{"E", LC(ExponentialCrossoverManager)},
		{"A", LC(ArithmeticCrossoverManager)},
});

std::vector<double> BinomialCrossoverManager::singleCrossover(std::vector<double>const& target, 
		std::vector<double>const& donor, double const Cr) const{
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

std::vector<double> ExponentialCrossoverManager::singleCrossover(std::vector<double>const& target, 
		std::vector<double>const& donor, double const Cr) const{
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

std::vector<double> ArithmeticCrossoverManager::singleCrossover(std::vector<double>const& target, 
		std::vector<double>const& donor, double const Cr) const{
	std::vector<double> x = target;
	double const k = rng.randDouble(0,1);
	std::vector<double> subtraction(D);
	subtract(donor, target, subtraction);
	scale(subtraction, k);
	add(x, subtraction, x);
	return x;
}
