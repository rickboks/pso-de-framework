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
#include "psode.h"
#include <limits>
#include <iostream>
#include <algorithm> 

PSODE::PSODE(UpdateManagerType const updateManagerType, 
	Topology topologyManagerType, Synchronicity const synchronicity, MutationType const mutationType, 
	CrossoverType const crossoverType, SelectionType selection,DEAdaptationType adaptionType):
		HybridAlgorithm(updateManagerType, topologyManagerType, synchronicity, 
			mutationType, crossoverType, selection, adaptionType){
}

PSODE::PSODE(hybrid_config config):
		HybridAlgorithm(config){
}

std::vector<Particle*> PSODE::copyPopulation(std::vector<Particle*>const& particles){
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

void PSODE::reset(){
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

PSODE::~PSODE(){
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

void PSODE::run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger,
    		int const evalBudget, int popSize, std::map<int,double> particleUpdateParams){

	this->problem=problem;
	this->logger=logger;
	if (config.synchronicity == SYNCHRONOUS)
		runSynchronous(evalBudget, popSize, particleUpdateParams);
	else if (config.synchronicity == ASYNCHRONOUS)
		runAsynchronous(evalBudget, popSize, particleUpdateParams);
}

void PSODE::runSynchronous(int const evalBudget, int popSize, 
	std::map<int,double> particleUpdateParams){

	std::vector<Particle*> p0;
	std::vector<Particle*> p1;
	std::vector<Particle*> p2;
		
	topologyManager = TopologyManager::createTopologyManager(config.topology, particles);
	popSize = topologyManager->getClosestValidPopulationSize(popSize);	

	int const D = problem->IOHprofiler_get_number_of_variables();
	std::vector<double> smallest = problem->IOHprofiler_get_lowerbound();
	std::vector<double> largest = problem->IOHprofiler_get_upperbound();

	ParticleUpdateSettings settings(config.update, particleUpdateParams, 
				smallest, largest);

	for (int i = 0; i < popSize; i++)
		particles.push_back(new Particle(D, settings));

	for (Particle* const p : particles)
		p->randomize(settings.xMax, settings.xMin);
	for (Particle* const p : particles)
		p->updateGbest();

	topologyManager->initialize();
	mutationManager = MutationManager<Particle>::createMutationManager(config.mutation, D);
	crossoverManager = CrossoverManager<Particle>::createCrossoverManager(config.crossover, D);
	adaptationManager = DEAdaptationManager::createDEAdaptationManager(config.adaptation);
	selectionManager = SelectionManager::createSelectionManager(config.selection, D, adaptationManager);

	int iterations = 0;
	double bestFitness = std::numeric_limits<double>::max();
	int notImproved = 0;
	bool improved;

	std::vector<double> Fs(popSize);
	std::vector<double> Crs(popSize);	

	while (	problem->IOHprofiler_get_evaluations() < evalBudget &&
			!problem->IOHprofiler_hit_optimal()){
		
		improved = false;

		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);
		adaptationManager->reset();

		for (int i = 0; i < popSize; i++){
			double y = particles[i]->evaluate(problem,logger);

			if (y < bestFitness){
				improved = true;
				bestFitness = y;
			}
		}

		improved ? notImproved = 0 : notImproved++;

		p0 = copyPopulation(particles);
		
		for (int i = 0; i < popSize; i++)
			p0[i]->updatePbest();
		for (int i = 0; i < popSize; i++)
			p0[i]->updateGbest();
		for (int i = 0; i < popSize; i++)
			p0[i]->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/double(evalBudget));

		p1 = mutationManager->mutate(particles, Fs);
		p2 = crossoverManager->crossover(particles,p1, Crs);
			
		improved = false;
		for (unsigned int i = 0; i < p2.size(); i++){
			p0[i]->evaluate(problem,logger);
			p2[i]->evaluate(problem,logger);
			
			double minF = std::min(p0[i]->getFitness(), p2[i]->getFitness());	
			if (minF < bestFitness){
				improved = true;
				bestFitness = minF;
			}
		}

		selectionManager->selection(particles, p0, p2);
		adaptationManager->update();

		improved ? notImproved = 0 : notImproved++;

		for (int i = 0; i < popSize; i++) {
			delete p0[i]; 
			delete p1[i];
			delete p2[i];
		}

		iterations++;	
		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/double(evalBudget));	
	}

	reset();
}

void PSODE::runAsynchronous(int const evalBudget, 
	int popSize, std::map<int,double> particleUpdateParams){

	std::vector<Particle*> p0;
	std::vector<Particle*> p1;
	std::vector<Particle*> p2;
		
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
	selectionManager = SelectionManager::createSelectionManager(config.selection, D, adaptationManager);

	int iterations = 0;

	double bestFitness = std::numeric_limits<double>::max();
	int notImproved = 0;
	bool improved;

	std::vector<double> Fs(popSize);
	std::vector<double> Crs(popSize);

	while (	problem->IOHprofiler_get_evaluations() < evalBudget &&
			!problem->IOHprofiler_hit_optimal()){
		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);
		adaptationManager->reset();
		
		improved = false;
		
		for (int i = 0; i < popSize; i++){
			double y = particles[i]->evaluate(problem,logger);

			if (y < bestFitness){
				improved = true;
				bestFitness = y;
			}
		}

		improved ? notImproved = 0 : notImproved++;

		p0 = copyPopulation(particles);
		
		for (int i = 0; i < popSize; i++){
			p0[i]->updatePbest();
			p0[i]->updateGbest();			
			p0[i]->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/double(evalBudget));
		}		

		p1 = mutationManager->mutate(particles, Fs);
		p2 = crossoverManager->crossover(particles,p1, Crs);
			
		improved = false;
		for (unsigned int i = 0; i < p2.size(); i++){
			p0[i]->evaluate(problem,logger);
			p2[i]->evaluate(problem,logger);


			double minF = std::min(p0[i]->getFitness(), p2[i]->getFitness());	
			if (minF < bestFitness){
				improved = true;
				bestFitness = minF;
			}
		}

		selectionManager->selection(particles, p0, p2);
		adaptationManager->update();

		improved ? notImproved = 0 : notImproved++;

		for (int i = 0; i < popSize; i++) {
			delete p0[i]; 
			delete p1[i];
			delete p2[i];
		}

		iterations++;	
		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/double(evalBudget));	
	}

	reset();
}

std::string PSODE::getIdString() const{
	return InstanceNamer::getName(config);
}
