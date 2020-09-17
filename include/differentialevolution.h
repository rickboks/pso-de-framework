#pragma once
#include "mutationmanager.h"
#include "particle.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "deadaptationmanager.h"

template <typename T>
class IOHprofiler_problem;

class IOHprofiler_csv_logger;

struct DEConfig {
	DEConfig(std::string const mutation, std::string const crossover, std::string const adaptation, std::string const constraintHandler)
	: mutation(mutation), crossover(crossover), adaptation(adaptation), constraintHandler(constraintHandler){}

	std::string const mutation, crossover, adaptation, constraintHandler;
};

class DifferentialEvolution {
	private:
		DEConfig const config;
		std::vector<Particle*> genomes;
		MutationManager* mutationManager;
		CrossoverManager* crossoverManager;
		DEAdaptationManager* adaptationManager;
		ConstraintHandler* deCH;
		int dimension;
		int popSize;
		std::shared_ptr<IOHprofiler_csv_logger> logger;
		std::shared_ptr<IOHprofiler_problem<double> > problem;
	public:
		DifferentialEvolution(DEConfig config);
		void run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger,
    		int const evalBudget, int const popSize);
		std::string getIdString() const;
};
