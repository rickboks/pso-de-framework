#pragma once
#include <random>
#include <stdexcept>
#include <vector>
#include<algorithm>
#include "rng.h"
#include "particle.h"

class CrossoverManager {
	protected:
		int const D;
	public:
		CrossoverManager(int const D);
		virtual ~CrossoverManager();

		std::vector<Solution*> crossover(std::vector<Solution*>const& genomes, std::vector<Solution*>const& mutants, std::vector<double>const& Crs) const;

		virtual std::vector<double> singleCrossover(std::vector<double>const& target, 
			std::vector<double>const& donor, double const Cr) const =0;
};

extern std::map<std::string, std::function<CrossoverManager* (int const)>> const crossovers;

class BinomialCrossoverManager : public CrossoverManager {
	public:
		BinomialCrossoverManager(int const D): CrossoverManager(D){};
		std::vector<double> singleCrossover(std::vector<double>const& target, std::vector<double>const& donor, double const Cr) const;
};

class ExponentialCrossoverManager : public CrossoverManager {
	public:
		ExponentialCrossoverManager(int const D): CrossoverManager(D){};
		std::vector<double> singleCrossover(std::vector<double>const& target, std::vector<double>const& donor, double const Cr) const;
};

class ArithmeticCrossoverManager : public CrossoverManager {
	public:
		ArithmeticCrossoverManager(int const D): CrossoverManager(D){};
		std::vector<double> singleCrossover(std::vector<double>const& target, std::vector<double>const& donor, double const Cr) const;
};
