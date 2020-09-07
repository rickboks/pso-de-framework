#include <IOHprofiler_problem.h>
#include <IOHprofiler_csv_logger.h>
#include "differentialevolution.h"
#include "rng.h"
#include "utilities.h"
#include "deadaptationmanager.h"
#include "instancenamer.h"
#include "util.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <functional>
#include <string>

DifferentialEvolution::DifferentialEvolution( MutationType const mutationType, CrossoverType const crossoverType, 
	DEAdaptationType const adaptationType)
	: mutationType(mutationType), crossoverType(crossoverType), adaptationType(adaptationType){
}

void DifferentialEvolution::run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger, 
    		int const evalBudget, int const popSize){

	this->logger = logger;
	this->problem = problem;
	this->popSize = popSize;

	dimension = problem->IOHprofiler_get_number_of_variables();

	for (int i = 0; i < popSize; i++)
		genomes.push_back(new Particle(dimension));

	for (Particle* const p : genomes)
		p->randomize(problem->IOHprofiler_get_lowerbound(), problem->IOHprofiler_get_upperbound());

	std::vector<Particle*> donors;
	std::vector<Particle*> trials;

	crossoverManager = CrossoverManager::createCrossoverManager(crossoverType, dimension);
	mutationManager = MutationManager::createMutationManager(mutationType, dimension);
	adaptationManager = DEAdaptationManager::createDEAdaptationManager(adaptationType);

	std::vector<double> Fs(popSize);
	std::vector<double> Crs(popSize);

	while (
		//noImprovement < 100 && 
		problem->IOHprofiler_get_evaluations() < evalBudget && 
		!problem->IOHprofiler_hit_optimal()){

		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);
		adaptationManager->reset();
		
		donors = mutationManager->mutate(genomes,Fs);
		trials = crossoverManager->crossover(genomes, donors, Crs);

		for (Particle* m : donors)
			delete m;
		
		for (int i = 0; i < popSize; i++){
			double parentF = genomes[i]->evaluate(problem,logger);
			double trialF = trials[i]->evaluate(problem,logger);

			if (trialF < parentF){
				genomes[i]->setX(trials[i]->getX(), trials[i]->getFitness(), false);
				adaptationManager->successfulIndex(i);				
			}
		}
		
		for (Particle* g : trials)
				delete g;

		adaptationManager->update();
	}

	for (Particle* d : genomes)
		delete d;
	
	delete mutationManager;
	delete crossoverManager;
	delete adaptationManager;

	genomes.clear();
}

std::string DifferentialEvolution::getIdString() const {
	return InstanceNamer::getName(mutationType, crossoverType, adaptationType);
}
