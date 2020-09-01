#include <IOHprofiler_problem.h>
#include <IOHprofiler_csv_logger.h>
#include "differentialevolution.h"
#include "deinitializer.h"
#include "rng.h"
#include "utilities.h"
#include "deadaptationmanager.h"
#include "instancenamer.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <functional>
#include <string>

bool comparePtrs(Genome* a, Genome* b);

DifferentialEvolution::DifferentialEvolution(DEInitializationType initializationType, 
	MutationType const mutationType, CrossoverType const crossoverType, 
	DEAdaptationType const adaptationType, bool const jumpOpposition)
	: mutationType(mutationType), crossoverType(crossoverType), 
	initializationType(initializationType), adaptationType(adaptationType), jumpOpposition(jumpOpposition){
}

void DifferentialEvolution::run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger, 
    		int const evalBudget, int const popSize){

	this->logger = logger;
	this->problem = problem;

	this->popSize = popSize;
	dimension = problem->IOHprofiler_get_number_of_variables();
	for (int i = 0; i < popSize; i++)
		genomes.push_back(new Genome(dimension));

	DEInitializer initializer(initializationType, problem, logger);
	initializer.initialize(genomes);
	double bestFitness = std::numeric_limits<double>::max();

	for (Genome* genome : genomes)
		if (genome->getFitness() < bestFitness)
			bestFitness = genome->getFitness();

	int noImprovement = 0;
	bool improved = false;

	std::vector<Genome*> donors;
	std::vector<Genome*> trials;

	crossoverManager = CrossoverManager<Genome>::createCrossoverManager(crossoverType, dimension);
	mutationManager = MutationManager<Genome>::createMutationManager(mutationType, dimension);
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
		improved = false;
		
		donors = mutationManager->mutate(genomes,Fs);
		trials = crossoverManager->crossover(genomes, donors, Crs);

		for (Genome* m : donors)
			delete m;
		
		for (int i = 0; i < popSize; i++){
			double parentF = genomes[i]->evaluate(problem,logger);
			if (parentF < bestFitness){
				improved = true;
				bestFitness = parentF;
			}

			double trialF = trials[i]->evaluate(problem,logger);
			if (trialF < bestFitness){
				improved = true;
				bestFitness = trialF;
			}

			if (trialF < parentF){
				genomes[i]->setPosition(trials[i]->getPosition(), trials[i]->getFitness());
				adaptationManager->successfulIndex(i);				
			}
		}
		
		for (Genome* g : trials){
				delete g;
		}

		if (jumpOpposition)
			oppositionGenerationJump();

		adaptationManager->update();
		improved ? noImprovement=0 : noImprovement++;
	}

	for (Genome* d : genomes)
		delete d;
	
	delete mutationManager;
	delete crossoverManager;
	delete adaptationManager;

	genomes.clear();
}

std::string DifferentialEvolution::getIdString() const {
	return InstanceNamer::getName(initializationType, mutationType, crossoverType, adaptationType, jumpOpposition);
}

void DifferentialEvolution::oppositionGenerationJump(){

	if (rng.randDouble(0,1) > 0.3){
		return;
	}

	std::vector<double> minValues(dimension, std::numeric_limits<double>::max());
	std::vector<double> maxValues(dimension, -std::numeric_limits<double>::max());

	for (Genome* genome: genomes){
		std::vector<double> x = genome->getPosition();

		for (int i = 0; i < dimension; i++){
			if (x[i] < minValues[i])
				minValues[i] = x[i];
			if (x[i] > maxValues[i])
				maxValues[i] = x[i];
		}
	}

	std::vector<Genome*> combined;
	combined.reserve(popSize * 2);

	for (Genome* g : genomes) {
		std::vector<double> x = g->getPosition();
		for (int i = 0; i < dimension; i++){
			x[i] = minValues[i] + maxValues[i] - x[i];
		}

		Genome* genome = new Genome(x);
		genome->evaluate(problem,logger);

		combined.push_back(genome);
	}

	combined.insert(combined.end(), genomes.begin(), genomes.end());

	std::sort(combined.begin(), combined.end(), comparePtrs);

	for (int i = 0; i < popSize; i++){
		genomes[i] = combined[i];
	}

	for (int i = popSize; i < (int)combined.size(); i++){
		delete combined[i];
	}
}
