#pragma once
#include "mutationmanager.h"
#include "particle.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "deadaptationmanager.h"

template <typename T>
class IOHprofiler_problem;

class IOHprofiler_csv_logger;

class DifferentialEvolution {
	private:
		std::vector<Particle*> genomes;
		MutationType const mutationType;
		CrossoverType const crossoverType;
		DEAdaptationType const adaptationType;
		MutationManager* mutationManager;
		CrossoverManager* crossoverManager;
		DEAdaptationManager* adaptationManager;
		int dimension;
		int popSize;
		std::shared_ptr<IOHprofiler_csv_logger> logger;
		std::shared_ptr<IOHprofiler_problem<double> > problem;
	public:
		DifferentialEvolution(MutationType const mutationType, CrossoverType const crossoverType, DEAdaptationType adaptationType);
		void run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger,
    		int const evalBudget, int const popSize);
		std::string getIdString() const;
};
