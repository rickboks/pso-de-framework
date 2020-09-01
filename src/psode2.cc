#include <IOHprofiler_problem.h>
#include <IOHprofiler_csv_logger.h>
#include "hybridalgorithm.h"
#include "particle.h"
#include "topologymanager.h"
#include "genome.h"
#include "particleupdatesettings.h"
#include "vectoroperations.h"
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

std::vector<Particle*> PSODE2::copyPopulation(std::vector<Particle*>const& particles){
	std::vector<Particle*> copy;
	copy.reserve(particles.size());

	std::map<Particle*, Particle*> newAddresses;	

	for (int i = 0; i < (int)particles.size(); i++){
		copy.push_back(new Particle(*(particles[i])));
		newAddresses[particles[i]] = copy[i];
	}

	for (Particle* p: copy){
		p->replaceNeighbors(newAddresses);
	}

	return copy;
}

void PSODE2::reset(){
	delete topologyManager;
	delete mutationManager;
	delete crossoverManager;
	delete adaptationManager;
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
	topologyManager = TopologyManager::createTopologyManager(config.topology, particles);
	popSize = topologyManager->getClosestValidPopulationSize(popSize);	

	int const D = problem->IOHprofiler_get_number_of_variables(); /// dimension
	std::vector<double> smallest = problem->IOHprofiler_get_lowerbound();
	std::vector<double> largest = problem->IOHprofiler_get_upperbound();

	ParticleUpdateSettings settings(config.update, particleUpdateParams, smallest, largest);

	for (int i = 0; i < popSize; i++)
		particles.push_back(new Particle(D, settings));

	for (Particle* const p : particles)
		p->randomize(settings.xMax, settings.xMin);

	topologyManager->initialize();
	mutationManager = MutationManager<Particle>::createMutationManager(config.mutation, D);
	crossoverManager = CrossoverManager<Particle>::createCrossoverManager(config.crossover, D);
	adaptationManager = DEAdaptationManager::createDEAdaptationManager(config.adaptation);
	//selectionManager = SelectionManager::createSelectionManager(config.selection, D, adaptationManager);

	int iterations = 0;
	double bestFitness = std::numeric_limits<double>::max();
	int notImproved = 0;
	bool improved;

	int split = particles.size() / 2;
	std::vector<int> deIndices(split);
	for (int i = 0; i < split; i++) 
		deIndices.push_back(i);
	std::vector<Particle*> psoPop = std::vector<Particle*>(particles.begin(), particles.begin() + split);

	std::vector<Particle*> dePop = std::vector<Particle*>(particles.begin() + split, particles.end());

	//Evaluate DE population to avoid uninitialized values in *best* mutation schemes
	for (Particle* p : dePop)
		p->evaluate(problem, logger);

	std::vector<double> Fs(dePop.size());
	std::vector<double> Crs(dePop.size());

	while (problem->IOHprofiler_get_evaluations() < evalBudget &&
			!problem->IOHprofiler_hit_optimal()){

		// Get new DE parameters from the adaptation manager (JADE or constant)
		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);
		adaptationManager->reset();

		// Measure fitness of the PSO population
		// Update pbest and gbest of PSO population
		// Update velocity and position of PSO population
		improved = false;
		for (int i = 0; i < psoPop.size(); i++){
			double y = particles[i]->evaluate(problem,logger);
			if (y < bestFitness){
				improved = true;
				bestFitness = y;
			}

			psoPop[i]->updatePbest();
			psoPop[i]->updateGbest();
			psoPop[i]->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/double(evalBudget));			
		}

		// Perform mutation and crossover
		std::vector<Particle*> donors = mutationManager->mutate(particles, deIndices, Fs);
		std::vector<Particle*> trials = crossoverManager->crossover(dePop, donors, Crs);

		// We don't need the donor vectors resulting from mutation anymore
		for (Particle* d : donors) 
			delete d;

		improved = false;
		for (int i = 0; i < dePop.size(); i++){
			//Evaluate the parent vector
			double parentF = dePop[i]->evaluate(problem,logger);
			if (parentF < bestFitness){
				improved = true;
				bestFitness = parentF;
			}

			//Evaluate the trial vector
			double trialF = trials[i]->evaluate(problem,logger);
			if (trialF < bestFitness){
				improved = true;
				bestFitness = trialF;
			}

			// Perform selection
			if (trials[i]->getFitness() < dePop[i]->getFitness()){
				dePop[i]->setPosition(trials[i]->getPosition(), trials[i]->getFitness());
				adaptationManager->successfulIndex(i);
			}

			dePop[i]->updatePbest();
			dePop[i]->updateGbest();
		}

		for (Particle* p : trials)
			delete p;

		adaptationManager->update();
		improved ? notImproved = 0 : notImproved++;
		iterations++;	
		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/evalBudget);	
	}
	reset();
}

std::string PSODE2::getIdString() const{
	return InstanceNamer::getName(config);
}
