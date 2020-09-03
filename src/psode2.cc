#include <IOHprofiler_problem.h>
#include <IOHprofiler_csv_logger.h>
#include "hybridalgorithm.h"
#include "repairhandler.h"
#include "particle.h"
#include "topologymanager.h"
#include "genome.h"
#include "particleupdatesettings.h"
#include "util.h"
#include "mutationmanager.h"
#include "instancenamer.h"
#include "psode2.h"
#include <limits>
#include <iostream>
#include <algorithm> 

PSODE2::PSODE2(UpdateManagerType const updateManagerType, 
	Topology topologyManagerType, Synchronicity const synchronicity, MutationType const mutationType, 
	CrossoverType const crossoverType, SelectionType selection,DEAdaptationType adaptionType):
		HybridAlgorithm(updateManagerType, topologyManagerType, synchronicity, 
			mutationType, crossoverType, selection, adaptionType){
}

PSODE2::PSODE2(hybrid_config config):
		HybridAlgorithm(config){
}

void PSODE2::reset(){
	delete topologyManager;
	delete mutationManager;
	delete crossoverManager;
	delete adaptationManager;
	delete constraintHandler;
	topologyManager = NULL;
	mutationManager = NULL;
	crossoverManager = NULL;
	adaptationManager = NULL;

	
	for (Particle* const particle : particles)
		delete particle;

	particles.clear();
}

PSODE2::~PSODE2(){
	if (topologyManager != NULL)
		delete topologyManager;
	if (mutationManager != NULL)
		delete mutationManager;
	if (crossoverManager != NULL)
		delete crossoverManager;
	if (adaptationManager != NULL)
		delete adaptationManager;

	if (!particles.empty())
		for (Particle* const particle : particles)
			delete particle;
}

void PSODE2::run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger,
    		int const evalBudget, int popSize, std::map<int,double> particleUpdateParams){

	this->problem=problem;
	this->logger=logger;
	//if (config.synchronicity == SYNCHRONOUS)
		//runSynchronous(evalBudget, popSize, particleUpdateParams);
	//if (config.synchronicity == ASYNCHRONOUS)
		runAsynchronous(evalBudget, popSize, particleUpdateParams);
}

void PSODE2::runAsynchronous(int const evalBudget, int popSize, std::map<int,double> particleUpdateParams){
	int const D = problem->IOHprofiler_get_number_of_variables(); /// dimension
	logging=false;

	std::vector<double> smallest = problem->IOHprofiler_get_lowerbound();
	std::vector<double> largest = problem->IOHprofiler_get_upperbound();

	ParticleUpdateSettings settings(config.update, particleUpdateParams, smallest, largest);

	for (int i = 0; i < popSize; i++)
		particles.push_back(new Particle(D, settings));

	for (Particle* const p : particles)
		p->randomize(settings.xMax, settings.xMin);

	mutationManager = MutationManager<Particle>::createMutationManager(config.mutation, D);
	crossoverManager = CrossoverManager<Particle>::createCrossoverManager(config.crossover, D);
	adaptationManager = DEAdaptationManager::createDEAdaptationManager(config.adaptation);
	constraintHandler = new ProjectionRepair(smallest, largest);

	int split = particles.size() / 2;
	psoPop = std::vector<Particle*>(particles.begin(), particles.begin() + split);
	dePop = std::vector<Particle*>(particles.begin() + split, particles.end());

	topologyManager = TopologyManager::createTopologyManager(config.topology, psoPop);
	topologyManager->initialize();

	for (Particle* p : particles)
		p->evaluate(problem, logger);

	std::vector<double> Fs(dePop.size());
	std::vector<double> Crs(dePop.size());

	logStart();
	int iterations = 0;
	while (problem->IOHprofiler_get_evaluations() < evalBudget &&
			!problem->IOHprofiler_hit_optimal()){

		// Get new DE parameters from the adaptation manager (JADE or constant)
		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);
		adaptationManager->reset();

		// Measure fitness of the PSO population
		// Update pbest and gbest of PSO population
		// Update velocity and position of PSO population
		for (int i = 0; i < psoPop.size(); i++){
			psoPop[i]->updatePbest();
			psoPop[i]->updateGbest();
			psoPop[i]->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/double(evalBudget));			
			psoPop[i]->repair(constraintHandler);
			psoPop[i]->evaluate(problem,logger);
		}

		// Perform mutation 
		std::vector<Particle*> donors = mutationManager->mutate(dePop, Fs);

		for (Particle* p : donors) 
			p->repair(constraintHandler); // Repair the mutants

		// Perform crossover
		std::vector<Particle*> trials = crossoverManager->crossover(dePop, donors, Crs);

		// We don't need the donor vectors resulting from mutation anymore
		for (Particle* d : donors) 
			delete d;

		for (int i = 0; i < dePop.size(); i++){
			//Evaluate the parent vector
			double parentF = dePop[i]->evaluate(problem,logger);

			//Evaluate the trial vector
			double trialF = trials[i]->evaluate(problem,logger);

			// Perform selection
			if (trialF < parentF){
				dePop[i]->setPosition(trials[i]->getPosition(), trials[i]->getFitness(), false);
				adaptationManager->successfulIndex(i);
			}
		}

		for (Particle* p : trials)
			delete p;

		if (iterations % 10 == 0)
			share();

		adaptationManager->update();
		iterations++;	
		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/evalBudget);	
		logPositions();
	}

	reset();
	logEnd();
}

std::string PSODE2::getIdString() const{
	return InstanceNamer::getName(config);
}

void PSODE2::logStart(){
	if (logging)
		std::cout << "LOGGER: START" << std::endl;
	// TODO: information about run
}

void PSODE2::logEnd(){
	if (logging)
		std::cout << "LOGGER: END" << std::endl;
}

// If PSO population has the best solution, replace the worst solution in DE with the best in PSO
// If DE population has the best solution, replace the worst solution in PSO with the best in DE
void PSODE2::share(){
	Particle* best_de = getBest(dePop);
	Particle* best_pso = getBest(psoPop);

	if (best_de->getFitness() < best_pso->getFitness()){
		Particle* worst_pso = getWorst(psoPop);
		worst_pso->setPosition(best_de->getPosition(), best_de->getFitness(), true); //Updates the velocity aswell
	} else {
		Particle* worst_de = getWorst(dePop);
		worst_de->setPosition(best_pso->getPosition(), best_pso->getFitness(), false); //Not here
	}
}

//void PSODE2::share(){
	//Particle* best_de = getBest(dePop);
	//Particle* best_pso = getBest(psoPop);

	//std::vector<double> x = best_pso->getPosition();
	//double y = best_pso->getFitness();
	//best_pso->setPosition(best_de->getPosition(), best_de->getFitness(), true); //Updates the velocity aswell
	//best_de->setPosition(x, y, false); //Not here
//}

void PSODE2::logPositions(){
	if (logging){
		std::cout << "LOGGER: " << std::endl;
		for (int i = 0; i < particles.size(); i++){
			std::cout << "LOGGER: ";
			std::vector<double> position = particles[i]->getPosition();
			for (int j = 0; j < position.size()-1 ; j++){
				std::cout << position[j] << " ";
			}
			std::cout << position[position.size()-1] << std::endl;
		}
	}
}
