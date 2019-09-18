#pragma once
#include <random>
#include <vector>

class Genome;
class CrossoverManager;

enum CrossoverType {
	BINOMIAL,
	EXPONENTIAL,
	CROSS_END
};

class CrossoverManagerFactory {
	public:
		static CrossoverManager* createCrossoverManager(CrossoverType const crossoverType, 
	std::vector<Genome*>& genomes, std::vector<Genome*>& mutants,double const Cr, int const D);
};

class CrossoverManager {
	protected:
		std::vector<Genome*>& genomes;
		std::vector<Genome*>& mutants;
		double const Cr;
		int const D;
	public:
		CrossoverManager(std::vector<Genome*>& genomes, std::vector<Genome*>& mutants, double const Cr, int const D);
		virtual ~CrossoverManager();
		virtual std::vector<Genome*> crossover() = 0;
};

class BinomialCrossoverManager : public CrossoverManager {
	private:

	public:
		BinomialCrossoverManager(std::vector<Genome*>& genomes, std::vector<Genome*>& mutants, double const Cr, int const D);
		std::vector<Genome*> crossover();
};

class ExponentialCrossoverManager : public CrossoverManager {
	private:

	public:
		ExponentialCrossoverManager(std::vector<Genome*>& genomes, std::vector<Genome*>& mutants, double const Cr, int const D);
		std::vector<Genome*> crossover();
};