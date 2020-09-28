#include <IOHprofiler_problem.h>
#include <IOHprofiler_csv_logger.h>
#include <random>
#include <algorithm>
#include <iostream>
#include <functional>
#include <string>
#include "differentialevolution.h"
#include "rng.h"
#include "utilities.h"
#include "deadaptationmanager.h"
#include "util.h"
#include "repairhandler.h"

DifferentialEvolution::DifferentialEvolution(DEConfig const config)
	: config(config){
}

void DifferentialEvolution::run(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger, 
    		int const evalBudget, int const popSize) const{

	int const D = problem->IOHprofiler_get_number_of_variables();
	std::vector<double> const lowerBound = problem->IOHprofiler_get_lowerbound();
	std::vector<double> const upperBound = problem->IOHprofiler_get_upperbound();

	std::vector<Solution*> genomes(popSize);
	for (int i = 0; i < popSize; i++){
		genomes[i] = new Solution(D);
		genomes[i]->randomize(lowerBound, upperBound);
		genomes[i]->evaluate(problem, logger);
	}

	DEConstraintHandler const* const deCH = deCHs.at(config.constraintHandler)(lowerBound, upperBound);
	CrossoverManager const* const crossoverManager = crossovers.at(config.crossover)(D);
	MutationManager* const mutationManager = mutations.at(config.mutation)(D, deCH);
	DEAdaptationManager* const adaptationManager = deAdaptations.at(config.adaptation)(popSize);

	std::vector<double> Fs(popSize);
	std::vector<double> Crs(popSize);

	while (problem->IOHprofiler_get_evaluations() < evalBudget && !problem->IOHprofiler_hit_optimal()){
		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);
		
		std::vector<Solution*> const donors = mutationManager->mutate(genomes,Fs);
		std::vector<Solution*> const trials = crossoverManager->crossover(genomes, donors, Crs);

		for (Solution* m : donors)
			delete m;
		
		std::vector<double> parentF(popSize), trialF(popSize);
		for (int i = 0; i < popSize; i++){
			parentF[i] = genomes[i]->getFitness();
			trialF[i] = trials[i]->evaluate(problem, logger);
			if (trialF[i] < parentF[i])
				genomes[i]->setX(trials[i]->getX(), trialF[i]);
		}
		
		for (Solution* g : trials)
			delete g;

		adaptationManager->update(parentF, trialF);
	}

	for (Solution* d : genomes)
		delete d;
	
	delete mutationManager;
	delete crossoverManager;
	delete adaptationManager;
	delete deCH;

	genomes.clear();
}

std::string DifferentialEvolution::getIdString() const {
	return "D_" + config.mutation + "_" + config.crossover + "_" + config.adaptation + "_" + config.constraintHandler;
}
