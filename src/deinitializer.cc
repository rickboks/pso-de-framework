#include "deinitializer.h"
#include "genome.h"
#include <random>
#include <vector>
#include <algorithm>
#include "iohsrc/Template/IOHprofiler_problem.hpp"

DEInitializer::DEInitializer(DEInitializationType const initializationType, std::shared_ptr<IOHprofiler_problem<double> > problem,
	std::shared_ptr<IOHprofiler_csv_logger> logger)
	: initializationType(initializationType), problem(problem), logger(logger){
}

void DEInitializer::initialize(std::vector<Genome*>& genomes) const{
	switch(initializationType){
		case RANDOM:
			initRandom(genomes);
			break;
		case OPPOSITION:
			initOpposition(genomes);
			break;
		default:
			throw std::invalid_argument("Error: Invalid DE Initialization type");
			break;
	}
}

void DEInitializer::initRandom(std::vector<Genome*>& genomes) const{
	for (Genome* g : genomes){
		g->randomize(problem->IOHprofiler_get_lowerbound(), problem->IOHprofiler_get_upperbound());
		g->evaluate(problem,logger);		
	}
}

bool comparePtrs(Genome* a, Genome* b) { 
	return (*a < *b); 
}

void DEInitializer::initOpposition(std::vector<Genome*>& genomes) const {
	int const popSize = genomes.size();
	int const D = problem->IOHprofiler_get_number_of_variables();
	std::vector<Genome*> combined;
	combined.reserve(genomes.size() * 2);

	std::vector<Genome*> random = genomes;
	for (Genome* g : random){
		g->randomize(problem->IOHprofiler_get_lowerbound(), problem->IOHprofiler_get_upperbound());
		g->evaluate(problem,logger);
	}

	std::vector<Genome*> opposite;
	opposite.reserve(random.size());

	for (Genome* g : random) {
		std::vector<double> x = g->getPosition();
		for (int i = 0; i < D; i++){
			x[i] = problem->IOHprofiler_get_lowerbound()[i] + problem->IOHprofiler_get_upperbound()[i] - x[i];
		}

		Genome* genome = new Genome(x);
		genome->evaluate(problem,logger);
		opposite.push_back(genome);
	}

	combined.insert(combined.end(), random.begin(), random.end());
	combined.insert(combined.end(), opposite.begin(), opposite.end());
	std::sort(combined.begin(), combined.end(), comparePtrs);
	
	for (int i = 0; i < popSize; i++){
		genomes[i] = combined[i];
	}

	for (int i = popSize; i < (int)combined.size(); i++){
		delete combined[i];
	}
}