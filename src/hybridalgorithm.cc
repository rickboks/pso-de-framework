#include "hybridalgorithm.h"
#include "particle.h"
#include "topologymanager.h"
#include "genome.h"
#include "particleupdatesettings.h"
#include "vectoroperations.h"
#include "coco.h"
#include "mutationmanager.h"
#include "instancenamer.h"
#include <limits>
#include <iostream>
#include <algorithm>   
#include <problem.h>

constexpr double DOUBLE_MIN = std::numeric_limits<double>::min();
constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();

HybridAlgorithm::HybridAlgorithm(UpdateManagerType const updateManagerType, 
			Topology topologyManagerType, Synchronicity const synchronicity, MutationType const mutationType, 
			CrossoverType const crossoverType, SelectionType selection,DEAdaptationType adaptionType):

	updateManagerType(updateManagerType),topologyManagerType(topologyManagerType), topologyManager(NULL), 
	synchronicity(synchronicity), mutationType(mutationType), crossoverType(crossoverType), 
	selectionType(selection), adaptationType(adaptionType), mutationManager(NULL), crossoverManager(NULL),
	adaptationManager(NULL){
}

std::vector<Particle*> HybridAlgorithm::copyPopulation(std::vector<Particle*>& particles){
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

void HybridAlgorithm::reset(){
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

HybridAlgorithm::~HybridAlgorithm(){
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

void HybridAlgorithm::run(Problem const problem, int const evalBudget, int popSize, std::map<int,double> particleUpdateParams, 
	double const F, double const Cr){
	if (synchronicity == SYNCHRONOUS)
		runSynchronous(problem, evalBudget, popSize, particleUpdateParams, F, Cr);
	else if (synchronicity == ASYNCHRONOUS){
		runAsynchronous(problem, evalBudget, popSize, particleUpdateParams, F, Cr);
	}
}

void HybridAlgorithm::runSynchronous(Problem const problem, int const evalBudget, int popSize, 
	std::map<int,double> particleUpdateParams, double const F, double const Cr){

	std::vector<Particle*> p1;
	std::vector<Particle*> p2;
	std::vector<Particle*> p3;
		
	topologyManager = TopologyManagerFactory::createTopologyManager(topologyManagerType, particles);
	popSize = topologyManager->getClosestValidPopulationSize(popSize);	

	int const D = coco_problem_get_dimension(problem.PROBLEM);
	double const* smallest = coco_problem_get_smallest_values_of_interest(problem.PROBLEM);
	double const* largest = coco_problem_get_largest_values_of_interest(problem.PROBLEM);

	ParticleUpdateSettings settings(updateManagerType, particleUpdateParams, 
				std::vector<double>(smallest, smallest + D), 
				std::vector<double>(largest, largest + D));

	double maxResV = 0;
	for (double d : settings.vMax)
		maxResV += d*d;

	maxResV = sqrt(maxResV);

	for (int i = 0; i < popSize; i++)
		particles.push_back(new Particle(coco_problem_get_dimension(problem.PROBLEM), settings));

	for (Particle* const p : particles)
		p->randomize(settings.xMax, settings.xMin);
	for (Particle* const p : particles)
		p->updateGbest();

	topologyManager->initialize();
	mutationManager = MutationManagerFactory::createMutationManager<Particle>(mutationType, D);
	crossoverManager = CrossoverManagerFactory::createCrossoverManager<Particle>(crossoverType, D);
	adaptationManager = DEAdaptationManagerFactory::createDEAdaptationManager(adaptationType);
	selectionManager = SelectionManagerFactory::createSelectionManager(selectionType, D, adaptationManager);

	int iterations = 0;
	double bestFitness = DOUBLE_MAX;
	int notImproved = 0;
	bool improved;
	int evaluations = 0;

	std::vector<double> Fs(popSize);
	std::vector<double> Crs(popSize);	

	while (	
			notImproved < 100 && 
			evaluations <= evalBudget &&
			!coco_problem_final_target_hit(problem.PROBLEM)){
		
		improved = false;

		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);
		adaptationManager->reset();

		for (int i = 0; i < popSize; i++){
			double y = particles[i]->evaluate(problem.evalFunc);
			evaluations++;

			if (y < bestFitness){
				improved = true;
				bestFitness = y;
			}
		}

		improved ? notImproved = 0 : notImproved++;

		p1 = copyPopulation(particles);
		
		for (int i = 0; i < popSize; i++)
			p1[i]->updatePbest();
		for (int i = 0; i < popSize; i++)
			p1[i]->updateGbest();
		for (int i = 0; i < popSize; i++)
			p1[i]->updateVelocityAndPosition(double(evaluations)/double(evalBudget));

		p2 = mutationManager->mutate(particles, Fs);
		p3 = crossoverManager->crossover(particles,p2, Crs);
			
		improved = false;
		for (unsigned int i = 0; i < p3.size(); i++){
			p1[i]->evaluate(problem.evalFunc);
			p3[i]->evaluate(problem.evalFunc);
			evaluations+=2;

			double minF = std::min(p1[i]->getFitness(), p3[i]->getFitness());	
			if (minF < bestFitness){
				improved = true;
				bestFitness = minF;
			}
		}

		selectionManager->selection(particles, p1, p3);
		adaptationManager->update();

		improved ? notImproved = 0 : notImproved++;



		for (int i = 0; i < popSize; i++) {
			//delete p0[i];  
			delete p1[i]; 
			delete p2[i];
			delete p3[i];
		}

		iterations++;	
		topologyManager->update(double(evaluations)/double(evalBudget));	
	}

	reset();
}

void HybridAlgorithm::runAsynchronous(Problem const problem, int const evalBudget, 
	int popSize, std::map<int,double> particleUpdateParams,
	double const F, double const Cr){

	std::vector<Particle*> p1;
	std::vector<Particle*> p2;
	std::vector<Particle*> p3;
		
	topologyManager = TopologyManagerFactory::createTopologyManager(topologyManagerType, particles);
	popSize = topologyManager->getClosestValidPopulationSize(popSize);	

	int const D = coco_problem_get_dimension(problem.PROBLEM);
	double const* smallest = coco_problem_get_smallest_values_of_interest(problem.PROBLEM);
	double const* largest = coco_problem_get_largest_values_of_interest(problem.PROBLEM);

	ParticleUpdateSettings settings(updateManagerType, particleUpdateParams, 
				std::vector<double>(smallest, smallest + D), 
				std::vector<double>(largest, largest + D));

	double maxResV = 0;
	for (double d : settings.vMax)
		maxResV += d*d;

	maxResV = sqrt(maxResV);

	for (int i = 0; i < popSize; i++)
		particles.push_back(new Particle(coco_problem_get_dimension(problem.PROBLEM), settings));

	for (Particle* const p : particles)
		p->randomize(settings.xMax, settings.xMin);
	for (Particle* const p : particles)
		p->updateGbest();

	topologyManager->initialize();
	mutationManager = MutationManagerFactory::createMutationManager<Particle>(mutationType,D);
	crossoverManager = CrossoverManagerFactory::createCrossoverManager<Particle>(crossoverType, D);
	adaptationManager = DEAdaptationManagerFactory::createDEAdaptationManager(adaptationType);
	selectionManager = SelectionManagerFactory::createSelectionManager(selectionType, D, adaptationManager);

	int iterations = 0;

	double bestFitness = DOUBLE_MAX;
	int notImproved = 0;
	bool improved;
	int evaluations = 0;

	std::vector<double> Fs(popSize);
	std::vector<double> Crs(popSize);

	while (	notImproved < 100 &&
			evaluations <= evalBudget &&
			!coco_problem_final_target_hit(problem.PROBLEM)){

		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);
		adaptationManager->reset();
		
		improved = false;
		p1 = copyPopulation(particles);
		for (int i = 0; i < popSize; i++){
			double y = particles[i]->evaluate(problem.evalFunc);
			evaluations++;

			if (y < bestFitness){
				improved = true;
				bestFitness = y;
			}

			p1[i]->updatePbest();
			p1[i]->updateGbest();			
			p1[i]->updateVelocityAndPosition(double(evaluations)/double(evalBudget));
		}

		improved ? notImproved = 0 : notImproved++;

		p2 = mutationManager->mutate(particles, Fs);
		p3 = crossoverManager->crossover(particles,p2, Crs);
			
		improved = false;
		for (unsigned int i = 0; i < p3.size(); i++){
			p1[i]->evaluate(problem.evalFunc);
			p3[i]->evaluate(problem.evalFunc);
			evaluations+=2;

			double minF = std::min(p1[i]->getFitness(), p3[i]->getFitness());	
			if (minF < bestFitness){
				improved = true;
				bestFitness = minF;
			}
		}

		selectionManager->selection(particles, p1, p3);
		adaptationManager->update();

		improved ? notImproved = 0 : notImproved++;

		for (int i = 0; i < popSize; i++) {
			//delete p0[i];  
			delete p1[i]; 
			delete p2[i];
			delete p3[i];
		}

		iterations++;	
		topologyManager->update(double(evaluations)/double(evalBudget));	
	}

	reset();
}

std::string HybridAlgorithm::getIdString() const{
	//std::string id = "H_";
	return InstanceNamer::getName(updateManagerType, topologyManagerType, synchronicity, mutationType, crossoverType,
		selectionType, adaptationType);
}