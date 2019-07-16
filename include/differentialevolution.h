#pragma once
#include "genome.h"
#include "problem.h"
#include "mutationmanager.h"
#include "deinitializer.h"
class DifferentialEvolution {
	private:
		std::vector<Genome*> genomes;
		MutationType const mutationType;
		CrossoverType const crossoverType;
		InitializationType const initializationType;
		bool const jumpOpposition;
		MutationManager* mutationManager;
		CrossoverManager* crossoverManager;
		int dimension;
		int popSize;

		void oppositionGenerationJump(Problem const problem);

	public:
		DifferentialEvolution(InitializationType initializationType, MutationType const mutationType, CrossoverType const crossoverType, bool const jumpOpposition);
		void run(Problem const problem, int const evalBudget, int const popSize, double const F, double const Cr);
		std::string getIdString() const;
};