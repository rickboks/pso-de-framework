#pragma once
#include "genome.h"
#include "problem.h"
#include "mutationmanager.h"
#include "deinitializer.h"
#include "deadaptationmanager.h"

class DifferentialEvolution {
	private:
		std::vector<Genome*> genomes;
		MutationType const mutationType;
		CrossoverType const crossoverType;
		DEInitializationType const initializationType;
		DEAdaptationType const adaptationType;
		bool const jumpOpposition;
		MutationManager<Genome>* mutationManager;
		CrossoverManager<Genome>* crossoverManager;
		DEAdaptationManager* adaptationManager;
		int dimension;
		int popSize;

		void oppositionGenerationJump(Problem const problem);

	public:
		DifferentialEvolution(DEInitializationType initializationType, MutationType const mutationType, 
			CrossoverType const crossoverType, DEAdaptationType adaptationType, bool const jumpOpposition);
		void run(Problem const problem, int const evalBudget, int const popSize, double const F, double const Cr);
		std::string getIdString() const;
};