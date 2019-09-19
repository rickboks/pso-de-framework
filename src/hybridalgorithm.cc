#include "hybridalgorithm.h"
#include "particle.h"
#include "topologymanager.h"
#include "genome.h"
#include "particleupdatesettings.h"
#include "vectoroperations.h"
#include "coco.h"
#include <limits>
#include <iostream>
#include <algorithm>   
#include <problem.h>

constexpr double DOUBLE_MIN = std::numeric_limits<double>::min();
constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();

void printVector(std::vector<double> v){
	std::cout << "{";
	for (unsigned int i = 0; i < v.size() -1; i++){
		std::cout << v[i] << " ";
	}

	std::cout << v[v.size() -1] << "}" << std::endl;
}

HybridAlgorithm::HybridAlgorithm(UpdateManagerType const updateManagerType, 
			Topology topologyManagerType, Synchronicity const synchronicity, MutationType const mutationType, 
			CrossoverType const crossoverType):

	updateManagerType(updateManagerType),topologyManagerType(topologyManagerType), topologyManager(NULL), 
	synchronicity(synchronicity), mutationType(mutationType), crossoverType(crossoverType), mutationManager(NULL), crossoverManager(NULL){
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
	topologyManager = NULL;
	mutationManager = NULL;
	crossoverManager = NULL;
	
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

void HybridAlgorithm::runSynchronous(Problem const problem, int const evalBudget, int popSize, std::map<int,double> particleUpdateParams,
	double const F, double const Cr){

	std::vector<Genome*> p0;
	std::vector<Particle*> p1;
	std::vector<Genome*> p2;
	std::vector<Genome*> p3;
		
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
	mutationManager = MutationManagerFactory::createMutationManager(mutationType, p0, F, D);
	
	int iterations = 0;

	double bestFitness = DOUBLE_MAX;
	int notImproved = 0;
	bool improved;
	int evaluations = 0;

	while (	
			notImproved < 100 && 
			evaluations <= evalBudget &&
			!coco_problem_final_target_hit(problem.PROBLEM)){
		
		improved = false;
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

		p0 = toGenomes(particles);
		p2 = mutationManager->mutate();	
			
		improved = false;
		for (unsigned int i = 0; i < p2.size(); i++){
			p1[i]->evaluate(problem.evalFunc);
			p2[i]->evaluate(problem.evalFunc);
			evaluations+=2;

			if (p1[i]->getFitness() < p2[i]->getFitness()){
				particles[i]->setPosition(p1[i]->getPosition(), p1[i]->getFitness());
				particles[i]->setVelocity(p1[i]->getVelocity());
			} else {
				std::vector<double> new_velocity(D);
				//reverse engineer velocity
				subtract(p2[i]->getPosition(), particles[i]->getPosition(), new_velocity);
				particles[i]->setVelocity(new_velocity);
				particles[i]->setPosition(p2[i]->getPosition(), p2[i]->getFitness());								
			}

			double minF = std::min(p1[i]->getFitness(), p2[i]->getFitness());	
			if (minF < bestFitness){
				improved = true;
				bestFitness = minF;
			}
		}

		improved ? notImproved = 0 : notImproved++;

		for (int i = 0; i < popSize; i++) {
			delete p0[i];  
			delete p1[i]; 
			delete p2[i];
		}

		iterations++;	
		topologyManager->update(double(evaluations)/double(evalBudget));	
	}

	reset();
}

void HybridAlgorithm::runAsynchronous(Problem const problem, int const evalBudget, int popSize, std::map<int,double> particleUpdateParams,
	double const F, double const Cr){

	std::vector<Genome*> p0;
	std::vector<Particle*> p1;
	std::vector<Genome*> p2;
	std::vector<Genome*> p3;
		
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
	mutationManager = MutationManagerFactory::createMutationManager(mutationType, p0, F, D);
	
	int iterations = 0;

	double bestFitness = DOUBLE_MAX;
	int notImproved = 0;
	bool improved;
	int evaluations = 0;

	while (	notImproved < 100 && 
			evaluations <= evalBudget &&
			!coco_problem_final_target_hit(problem.PROBLEM)){
		
		improved = false;
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

		p1 = copyPopulation(particles);
		


		p0 = toGenomes(particles);
		p2 = mutationManager->mutate();	
			
		improved = false;
		for (unsigned int i = 0; i < p2.size(); i++){
			p1[i]->evaluate(problem.evalFunc);
			p2[i]->evaluate(problem.evalFunc);
			evaluations+=2;

			if (p1[i]->getFitness() < p2[i]->getFitness()){
				particles[i]->setPosition(p1[i]->getPosition(), p1[i]->getFitness());
				particles[i]->setVelocity(p1[i]->getVelocity());
			} else {
				std::vector<double> new_velocity(D);
				//reverse engineer velocity
				subtract(p2[i]->getPosition(), particles[i]->getPosition(), new_velocity);
				particles[i]->setVelocity(new_velocity);
				particles[i]->setPosition(p2[i]->getPosition(), p2[i]->getFitness());								
			}

			double minF = std::min(p1[i]->getFitness(), p2[i]->getFitness());	
			if (minF < bestFitness){
				improved = true;
				bestFitness = minF;
			}
		}

		improved ? notImproved = 0 : notImproved++;

		for (int i = 0; i < popSize; i++) {
			delete p0[i];  
			delete p1[i]; 
			delete p2[i];
		}

		iterations++;	
		topologyManager->update(double(evaluations)/double(evalBudget));	
	}

	reset();
}

std::vector<Genome*> HybridAlgorithm::toGenomes(std::vector<Particle*>& particles){
	std::vector<Genome*> genomes;
	genomes.reserve(particles.size());
	for (Particle* p : particles){
		genomes.push_back(new Genome(p));
	}

	return genomes;
}
std::string HybridAlgorithm::getIdString() const{
	std::string id = "H_";

	switch (updateManagerType){
		case UpdateManagerType::INERTIA_WEIGHT:
		id += "I";
		break;
		case UpdateManagerType::DECR_INERTIA_WEIGHT:
		id += "D";
		break;
		case UpdateManagerType::VMAX:
		id += "V";
		break;
		case UpdateManagerType::CONSTRICTION_COEFFICIENT:
		id += "C";
		break;
		case UpdateManagerType::FIPS:
		id += "F";
		break;
		case UpdateManagerType::BARE_BONES:
		id += "B";
		break;
		default:
		id+="ERR";
		break;
	};

	id += "_";
	
	switch (topologyManagerType){
		case Topology::LBEST:
		id += "L";
		break;
		case Topology::GBEST:
		id += "G";
		break;
		case Topology::RANDOM_GRAPH:
		id += "R";
		break;
		case Topology::VON_NEUMANN:
		id += "N";
		break;
		case Topology::WHEEL:
		id += "W";
		break;
		case Topology::INCREASING:
		id += "I";
		break;
		case Topology::DECREASING:
		id += "D";
		break;
		case Topology::MULTI_SWARM:
		id += "M";
		break;
		default:
		id+="ERR";
		break;
	};

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

	return id;
}
