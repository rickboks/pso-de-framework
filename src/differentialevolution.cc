#include <IOHprofiler_problem.h>
#include <IOHprofiler_csv_logger.h>
#include <limits>
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
    		int const evalBudget, int const popSize) const {

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

	Logger logger("scratch/extra_data/" + getIdString() + ".dat");
	Logger loggerParams("scratch/extra_data/" + getIdString() + ".par");
	//Logger loggerAnimation("scratch/animations/" + getIdString() + "_f" +
			//std::to_string(problem->IOHprofiler_get_problem_id()) + "D" + std::to_string(D) + 
			//".log");

	loggerParams.start(problem->IOHprofiler_get_problem_id(), D);

	int iteration = 0;
	while (problem->IOHprofiler_get_evaluations() < evalBudget && !problem->IOHprofiler_hit_optimal()){
		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);

		if (iteration % 10 == 0)
			loggerParams.log(Fs, Crs);
		
		std::vector<Solution*> const donors = mutationManager->mutate(genomes,Fs);
		std::vector<Solution*> const trials = crossoverManager->crossover(genomes, donors, Crs);

		for (Solution* m : donors)
			delete m;

		std::vector<double> parentF(popSize), trialF(popSize);
		for (int i = 0; i < popSize; i++){
			parentF[i] = genomes[i]->getFitness();

			trials[i]->evaluate(problem, iohLogger);

			deCH->penalize(trials[i]); // This is done after and not before the evaluation, because otherwise it could loop endlessly

			trialF[i] = trials[i]->getFitness();

			int const numEval = problem->IOHprofiler_get_evaluations();
			if (numEval != 0 && numEval % 100000 == 0)
				percCorrected.push_back(double(deCH->getCorrections()) / numEval);

			if (trialF[i] < parentF[i])
				genomes[i]->setX(trials[i]->getX(), trialF[i]);
		}

		for (Solution* g : trials)
			delete g;

		//loggerAnimation.log(genomes);

		adaptationManager->update(parentF, trialF);
		iteration++;
	}

	if (percCorrected.empty()){
		double const perc = double(deCH->getCorrections()) / problem->IOHprofiler_get_evaluations();
		percCorrected.resize(3, perc);
	} else if (percCorrected.size() < 3){
		int const lastIndex = percCorrected.size() -1;
		for (int i = lastIndex+1; i < 3; i++)
			percCorrected.push_back(percCorrected[lastIndex]);
	}

	Solution const*const best = getBest(genomes);

	logger.log(problem->IOHprofiler_get_problem_id(), D, percCorrected, best->getX(), best->getFitness(), problem->IOHprofiler_get_evaluations());
	loggerParams.newLine();

	for (Solution* d : genomes)
		delete d;

	delete mutationManager;
	delete crossoverManager;
	delete adaptationManager;
	delete deCH;

	genomes.clear();
}

std::string DifferentialEvolution::getIdString() const {
	return /*"DE_" +*/ config.mutation + "_" + config.crossover + "_" /*+ config.adaptation + "_"*/ + config.constraintHandler;
}
