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

		virtual std::vector<Particle*> crossover(std::vector<Particle*>const& targets, 
			std::vector<Particle*>const& donors, std::vector<double>const& Crs) const = 0;

		virtual std::vector<double> singleCrossover(std::vector<double>const& target, 
			std::vector<double>const& donor, double const Cr) const =0;
};

extern std::map<std::string, std::function<CrossoverManager* (int const)>> const crossovers;

class BinomialCrossoverManager : public CrossoverManager {
	public:
		BinomialCrossoverManager(int const D);
		std::vector<Particle*> crossover(std::vector<Particle*>const& genomes, std::vector<Particle*>const& mutants, std::vector<double>const& Crs) const;
		std::vector<double> singleCrossover(std::vector<double>const& target, std::vector<double>const& donor, double const Cr) const;
};

class ExponentialCrossoverManager : public CrossoverManager {
	public:
		ExponentialCrossoverManager(int const D);
		std::vector<Particle*> crossover(std::vector<Particle*>const& genomes, std::vector<Particle*>const& mutants, std::vector<double>const& Crs) const;
		std::vector<double> singleCrossover(std::vector<double>const& target, std::vector<double>const& donor, double const Cr) const;
};
