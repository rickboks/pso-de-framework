#include "differentialevolution.h"
#include "deinitializer.h"
#include "coco.h"
#include "rng.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <functional>
#include <string>

bool comparePtrs(Genome* a, Genome* b);

DifferentialEvolution::DifferentialEvolution(InitializationType initializationType, MutationType const mutationType, 
	CrossoverType const crossoverType, bool const jumpOpposition)
	: mutationType(mutationType), crossoverType(crossoverType), 
	initializationType(initializationType), jumpOpposition(jumpOpposition){
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

	std::vector<Genome*> mutants;
	std::vector<Genome*> donors;
	std::vector<Genome*> newPopulation;
	std::vector<Genome*> oldPopulation;
	newPopulation.reserve(popSize);
	oldPopulation.reserve(popSize);

	crossoverManager = CrossoverManagerFactory::createCrossoverManager(crossoverType, genomes, mutants, Cr);
	mutationManager = MutationManagerFactory::createMutationManager(mutationType, genomes, F);

	while (noImprovement < 100 && evaluations < evalBudget && !coco_problem_final_target_hit(problem.PROBLEM)){
		improved = false;
		
		mutants = mutationManager->mutate();
		donors = crossoverManager->crossover();

		for (Genome* m : mutants)
			delete m;
		
		for (int i = 0; i < popSize; i++){
			double parentF = genomes[i]->evaluate(problem.evalFunc);
			if (parentF < bestFitness){
				improved = true;
				bestFitness = parentF;
			}
			evaluations++;
			if (evaluations >= evalBudget || coco_problem_final_target_hit(problem.PROBLEM)){
				for (int k = i; k < popSize; k++){
					oldPopulation.push_back(donors[k]);
					oldPopulation.push_back(genomes[k]);
				}
				goto stop;
			}

			double donorF = donors[i]->evaluate(problem.evalFunc);
			if (donorF < bestFitness){
				improved = true;
				bestFitness = donorF;
			}
			evaluations++;

			if (evaluations >= evalBudget || coco_problem_final_target_hit(problem.PROBLEM)){
				for (int k = i; k < popSize; k++){
					oldPopulation.push_back(donors[k]);
					oldPopulation.push_back(genomes[k]);
				}
				goto stop;
			}

			if (parentF < donorF){
				newPopulation.push_back(genomes[i]);
				oldPopulation.push_back(donors[i]);
			} else {
				newPopulation.push_back(donors[i]);
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

		improved ? noImprovement=0 : noImprovement++;		
	}

	for (Genome* d : genomes)
		delete d;
	
stop:
	delete mutationManager;
	delete crossoverManager;

	for (Genome* d : oldPopulation)
		delete d;

	for (Genome* d : newPopulation)
		delete d;

	

	genomes.clear();
}

std::string DifferentialEvolution::getIdString() const {
	std::string id = "DE_";

	switch(initializationType){
		case InitializationType::RANDOM:
			id += "R";
			break;
		case InitializationType::OPPOSITION:
			id += "O";
			break;
		default:
			id += "ERR";
			break;
	}

	id += "_"; 

	switch (mutationType){
		case MutationType::RAND_1:
			id += "R1";
			break;
		case MutationType::BEST_1:
			id += "B1";
			break;
		case MutationType::TTB_1:
			id += "TTB1";
			break;
		case MutationType::BEST_2:
			id += "B2";
			break;
		case MutationType::RAND_2:
			id += "R2";
			break;
		case MutationType::RAND_2_DIR:
			id += "R2D";
			break;
		case MutationType::NSDE:
			id += "NS";
			break;
		case MutationType::TOPOLOGY:
			id += "TOP";
			break;
		default:
			id += "ERR";
			break;
	}

	id += "_";

	switch (crossoverType){
		case CrossoverType::BINOMIAL:
			id += "BIN";
			break;
		case CrossoverType::EXPONENTIAL:
			id += "EXP";
			break;
		default:
			id += "ERR";
			break;
	}

	id += "_";

	switch (jumpOpposition){
		case true:
		 	id += "1";
		 	break;
		case false:
			id += "0";
			break;
	}

	return id;

}

void DifferentialEvolution::oppositionGenerationJump(Problem const problem){

	if (rng.randDouble(0,1) > 0.3){
		return;
	}

	std::vector<double> minValues(dimension, std::numeric_limits<double>::max());
	std::vector<double> maxValues(dimension, -std::numeric_limits<double>::max());

	for (Genome* genome: genomes){
		std::vector<double> x = genome->getX();

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
		std::vector<double> x = g->getX();
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