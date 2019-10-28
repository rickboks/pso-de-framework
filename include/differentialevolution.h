#pragma once
#include "genome.h"
#include "mutationmanager.h"
#include "deinitializer.h"
#include "deadaptationmanager.h"

template <typename T>
class IOHprofiler_problem;

class IOHprofiler_csv_logger;

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
		std::shared_ptr<IOHprofiler_csv_logger> logger;
		std::shared_ptr<IOHprofiler_problem<double> > problem;

		void oppositionGenerationJump();

	public:
		DifferentialEvolution(DEInitializationType initializationType, MutationType const mutationType, 
			CrossoverType const crossoverType, DEAdaptationType adaptationType, bool const jumpOpposition);
		void run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger,
    		int const evalBudget, int const popSize);
		std::string getIdString() const;
};