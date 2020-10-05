#include <IOHprofiler_problem.h>
#include <IOHprofiler_csv_logger.h>
#include <random>
#include <algorithm>
#include <iostream>
#include <functional>
#include <string>
#include <fstream> 
#include "differentialevolution.h"
#include "rng.h"
#include "utilities.h"
#include "deadaptationmanager.h"
#include "util.h"
#include "repairhandler.h"
#include "logger.h"

DifferentialEvolution::DifferentialEvolution(DEConfig const config)
	: config(config){
}

void DifferentialEvolution::run(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const iohLogger, 
    		int const evalBudget, int const popSize) const{

	int const D = problem->IOHprofiler_get_number_of_variables();
	std::vector<double> const lowerBound = problem->IOHprofiler_get_lowerbound();
	std::vector<double> const upperBound = problem->IOHprofiler_get_upperbound();

	std::vector<Solution*> genomes(popSize);
	for (int i = 0; i < popSize; i++){
		genomes[i] = new Solution(D);
		genomes[i]->randomize(lowerBound, upperBound);
		genomes[i]->evaluate(problem, iohLogger);
	}

	DEConstraintHandler * const deCH = deCHs.at(config.constraintHandler)(lowerBound, upperBound);
	CrossoverManager const* const crossoverManager = crossovers.at(config.crossover)(D);
	MutationManager* const mutationManager = mutations.at(config.mutation)(D, deCH);
	DEAdaptationManager* const adaptationManager = deAdaptations.at(config.adaptation)(popSize);

	std::vector<double> Fs(popSize);
	std::vector<double> Crs(popSize);
	std::vector<double> percCorrected; 

	std::filebuf fb;
	fb.open("test.txt", std::ios::out);
	Logger logger(&fb);

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
			trialF[i] = trials[i]->evaluate(problem, iohLogger);

			int const numEval = problem->IOHprofiler_get_evaluations();
			if (numEval != 0 && numEval % 100000 == 0)
				percCorrected.push_back(double(deCH->getCorrections()) / numEval); // save #corr every 100000 evals

			if (trialF[i] < parentF[i])
				genomes[i]->setX(trials[i]->getX(), trialF[i]);
		}

		logger.log(genomes);
		
		for (Solution* g : trials)
			delete g;

		adaptationManager->update(parentF, trialF);
	}

	//Solution const* const best = getBest(genomes);
	//std::cout << std::endl << getIdString() << "p" << popSize << "D" << D << "f" << problem->IOHprofiler_get_problem_id();
	//for (double d : percCorrected)
		//std::cout << d << " ";
	//for (double d : best->getX())
		//std::cout << d << " ";
	//std::cout << best->getFitness() << " " << problem->IOHprofiler_get_optimal()[0] << std::endl;

	for (Solution* d : genomes)
		delete d;
	delete mutationManager;
	delete crossoverManager;
	delete adaptationManager;
	delete deCH;

	genomes.clear();
}

std::string DifferentialEvolution::getIdString() const {
	return "DE_" + config.mutation + "_" + config.crossover + "_" + config.adaptation + "_" + config.constraintHandler;
}
