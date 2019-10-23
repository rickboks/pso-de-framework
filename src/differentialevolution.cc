#include "differentialevolution.h"
#include "deinitializer.h"
#include "coco.h"
#include "rng.h"
#include "utilities.h"
#include "deadaptationmanager.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <functional>
#include <string>

bool comparePtrs(Genome* a, Genome* b);

DifferentialEvolution::DifferentialEvolution(InitializationType initializationType, 
	MutationType const mutationType, CrossoverType const crossoverType, 
	DEAdaptationType const adaptationType, bool const jumpOpposition)
	: mutationType(mutationType), crossoverType(crossoverType), 
	initializationType(initializationType), adaptationType(adaptationType), jumpOpposition(jumpOpposition){
}

void DifferentialEvolution::run(Problem const problem, int const evalBudget, int const popSize, double const F, double const Cr){
	this->popSize = popSize;
	dimension = coco_problem_get_dimension(problem.PROBLEM);
	for (int i = 0; i < popSize; i++)
		genomes.push_back(new Genome(dimension));

	int evaluationsBefore = coco_problem_get_evaluations(problem.PROBLEM);
	DEInitializer initializer(initializationType, problem);
	initializer.initialize(genomes);
	int evaluations = coco_problem_get_evaluations(problem.PROBLEM) - evaluationsBefore;

	double bestFitness = std::numeric_limits<double>::max();

	for (Genome* genome : genomes)
		if (genome->getFitness() < bestFitness)
			bestFitness = genome->getFitness();

	int noImprovement = 0;
	bool improved = false;

	std::vector<Genome*> donors;
	std::vector<Genome*> trials;
	std::vector<Genome*> newPopulation;
	std::vector<Genome*> oldPopulation;
	newPopulation.reserve(popSize);
	oldPopulation.reserve(popSize);

	crossoverManager = CrossoverManagerFactory::createCrossoverManager<Genome>(crossoverType, dimension);
	mutationManager = MutationManagerFactory::createMutationManager<Genome>(mutationType, dimension);
	adaptationManager = DEAdaptationManagerFactory::createDEAdaptationManager(adaptationType);

	std::vector<double> Fs(popSize);
	std::vector<double> Crs(popSize);

	while (
		noImprovement < 100 && 
		evaluations < evalBudget && 
		!coco_problem_final_target_hit(problem.PROBLEM)){

		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);
		adaptationManager->reset();
		improved = false;
		
		donors = mutationManager->mutate(genomes,Fs);
		trials = crossoverManager->crossover(genomes, donors, Crs);

		for (Genome* m : donors)
			delete m;
		
		for (int i = 0; i < popSize; i++){
			double parentF = genomes[i]->evaluate(problem.evalFunc);
			if (parentF < bestFitness){
				improved = true;
				bestFitness = parentF;
			}
			evaluations++;

			double trialF = trials[i]->evaluate(problem.evalFunc);
			if (trialF < bestFitness){
				improved = true;
				bestFitness = trialF;
			}
			evaluations++;

			if (parentF < trialF){
				newPopulation.push_back(genomes[i]);
				oldPopulation.push_back(trials[i]);
			} else {
				adaptationManager->successfulIndex(i);				
				newPopulation.push_back(trials[i]);
				oldPopulation.push_back(genomes[i]);
			}
		}

		genomes = newPopulation;
		newPopulation.clear();

		for (Genome* g : oldPopulation)
			delete g;

		oldPopulation.clear();

		if (jumpOpposition)
			oppositionGenerationJump(problem);

		adaptationManager->update();
		improved ? noImprovement=0 : noImprovement++;
	}

	for (Genome* d : genomes)
		delete d;
	
	delete mutationManager;
	delete crossoverManager;
	delete adaptationManager;

	for (Genome* d : oldPopulation)
		delete d;

	for (Genome* d : newPopulation)
		delete d;	

	genomes.clear();
}

std::string DifferentialEvolution::getIdString() const {
	std::string id = "";

	switch(initializationType){
		case InitializationType::RANDOM: id += "R";	break;
		case InitializationType::OPPOSITION: id += "O";	break;
		default: id += "ERR"; break;
	}

	// id += "_"; 

	switch (mutationType){
		case MutationType::RAND_1: id += "R1"; break;
		case MutationType::BEST_1: id += "B1"; break;
		case MutationType::TTB_1: id += "T1"; break;
		case MutationType::BEST_2: id += "B2"; break;
		case MutationType::RAND_2: id += "R2"; break;
		case MutationType::RAND_2_DIR: id += "RD"; break;
		case MutationType::NSDE: id += "NS"; break;
		case MutationType::TOPOLOGY: id += "TOP"; break;
		case MutationType::TRIGONOMETRIC: id += "TR"; break;
		case MutationType::TTPB_1: id+= "PB"; break;
		case MutationType::TO1: id+= "O1"; break;
		case MutationType::TO2: id+= "O2"; break;
		default: id += "ERR"; break;
	}

	// id += "_";

	switch (crossoverType){
		case CrossoverType::BINOMIAL: id += "B"; break;
		case CrossoverType::EXPONENTIAL: id += "E"; break;
		default: id += "ERR"; break;
	}

	switch (adaptationType){
		case DEAdaptationType::JADE: id += "J"; break;
		case DEAdaptationType::NO: id += "N"; break;
		default: id += "ERR"; break;
	}

	// id += "_";

	jumpOpposition ? id += "1" : id += "0";

	return id;

}

void DifferentialEvolution::oppositionGenerationJump(Problem const problem){

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
		genome->evaluate(problem.evalFunc);
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