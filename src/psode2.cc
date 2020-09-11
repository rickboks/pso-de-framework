#include <IOHprofiler_problem.h>
#include <IOHprofiler_csv_logger.h>
#include "hybridalgorithm.h"
#include "repairhandler.h"
#include "particle.h"
#include "topologymanager.h"
#include "particleupdatesettings.h"
#include "util.h"
#include "mutationmanager.h"
#include "instancenamer.h"
#include "penaltyhandler.h"
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
	delete deCH;
	delete psoCH;

	for (Particle* const particle : particles)
		delete particle;

	particles.clear();
	dePop.clear();
	psoPop.clear();
}

PSODE2::~PSODE2(){}

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

	deCH = new ResamplingRepair(smallest, largest);
	psoCH = new DeathPenalty(smallest, largest);

	mutationManager = MutationManager::createMutationManager(config.mutation, D, deCH);
	crossoverManager = CrossoverManager::createCrossoverManager(config.crossover, D);
	adaptationManager = DEAdaptationManager::createDEAdaptationManager(config.adaptation);

	int split = popSize / 2;

	for (int i = 0; i < split; i++) psoPop.push_back(new Particle(D, settings));
	for (int i = split; i < popSize; i++) dePop.push_back(new Particle(D));
	
	// append the two populations
	particles = psoPop; 
	particles.insert(particles.end(), dePop.begin(), dePop.end());

	for (Particle* const p : particles)
		p->randomize(settings.xMax, settings.xMin);

	for (Particle* const p : particles)
		p->evaluate(problem, logger);

	topologyManager = TopologyManager::createTopologyManager(config.topology, psoPop);
	topologyManager->initialize();

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
		for (Particle* const p : psoPop){
			p->updatePbest();
			p->updateGbest();
			p->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/double(evalBudget));			
			psoCH->repair(p); // generic repair
			p->evaluate(problem,logger);
		}

		// Perform mutation 
		//
		std::vector<Particle*> donors = mutationManager->mutate(dePop, Fs);

		for (Particle* const p : donors) 
			psoCH->repair(p); // generic repair

		// Perform crossover
		std::vector<Particle*> trials = crossoverManager->crossover(dePop, donors, Crs);

		for (Particle* const d : donors) delete d;

		for (int i = 0; i < dePop.size(); i++){
			//Evaluate the parent vector
			dePop[i]->evaluate(problem,logger);

			//Evaluate the trial vector
			trials[i]->evaluate(problem,logger);

			// Perform selection
			if (trials[i]->getFitness() < dePop[i]->getFitness()){
				dePop[i]->setX(trials[i]->getX(), trials[i]->getFitness(), false);
				adaptationManager->successfulIndex(i);
			}
		}

		for (Particle* const p : trials)
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

void PSODE2::share(){
	Particle* best_de = getPBest(dePop, 0.1);
	Particle* best_pso = getPBest(psoPop, 0.1);

	std::vector<double> x = best_pso->getX();
	double y = best_pso->getFitness();

	best_pso->setX(best_de->getX(), best_de->getFitness(), true); //Updates the velocity as well
	best_de->setX(x, y, false); //Not here
}

void PSODE2::logPositions(){
	if (logging){
		std::cout << "LOGGER: " << std::endl;
		for (int i = 0; i < particles.size(); i++){
			std::cout << "LOGGER: ";
			std::vector<double> position = particles[i]->getX();
			for (int j = 0; j < position.size()-1 ; j++){
				std::cout << position[j] << " ";
			}
			std::cout << position[position.size()-1] << std::endl;
		}
	}
}
