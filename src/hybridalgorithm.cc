#include "hybridalgorithm.h"
#include "particle.h"
#include "topologymanager.h"
#include "genome.h"
#include "particleupdatesettings.h"
#include "coco.h"
#include <limits>
#include <iostream>
#include <problem.h>

constexpr double DOUBLE_MIN = std::numeric_limits<double>::min();
constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();

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
		throw std::invalid_argument("Error: Asynchronous updates not implemented in HybridAlgorithm");
	}

}

void HybridAlgorithm::runSynchronous(Problem const problem, int const evalBudget, int popSize, std::map<int,double> particleUpdateParams,
	double const F, double const Cr){

	std::vector<Genome*> p0;
	std::vector<Genome*> p1;
	std::vector<Genome*> p2;
	std::vector<Genome*> p3;
		
	topologyManager = TopologyManagerFactory::createTopologyManager(topologyManagerType, particles);
	popSize = topologyManager->getClosestValidPopulationSize(popSize);	

	int const dimension = coco_problem_get_dimension(problem.PROBLEM);
	double const* smallest = coco_problem_get_smallest_values_of_interest(problem.PROBLEM);
	double const* largest = coco_problem_get_largest_values_of_interest(problem.PROBLEM);

	ParticleUpdateSettings settings(updateManagerType, particleUpdateParams, 
				std::vector<double>(smallest, smallest + dimension), 
				std::vector<double>(largest, largest + dimension));

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

	p0 = toGenomes(particles);
	p1 = toGenomes(particles);

	crossoverManager = CrossoverManagerFactory::createCrossoverManager(crossoverType, p1, p2, Cr);
	mutationManager = MutationManagerFactory::createMutationManager(mutationType, p0, F);

	for (Genome* g : p0) delete g;
	for (Genome* g : p1) delete g;	
	
	int iterations = 0;

	double bestFitness = DOUBLE_MAX;
	int notImproved = 0;
	bool improved;
	int evaluations = 0;
	//Main loop
	while (	notImproved < 100 && 
			evaluations <= evalBudget &&
			!coco_problem_final_target_hit(problem.PROBLEM)){
		
		improved = false;

		for (int i = 0; i < popSize; i++){
			double y = particles[i]->evaluate(problem.evalFunc);
			evaluations++;
			if (evaluations >= evalBudget || coco_problem_final_target_hit(problem.PROBLEM)){
				reset();
				return;
			}

			if (y < bestFitness){
				improved = true;
				bestFitness = y;
			}
		}

		improved ? notImproved = 0 : notImproved++;

		std::vector<Particle*> copy = copyPopulation(particles);
		
		for (int i = 0; i < popSize; i++)
			copy[i]->updatePbest();

		for (int i = 0; i < popSize; i++)
			copy[i]->updateGbest();

		for (int i = 0; i < popSize; i++)
			copy[i]->updateVelocityAndPosition(double(evaluations)/double(evalBudget));

		p0 = toGenomes(particles);
		p1 = toGenomes(copy);
		p2 = mutationManager->mutate();
		p3 = crossoverManager->crossover();

		for (Genome* g : p1) delete g;
		for (Genome* g : p2) delete g;
		for (Particle* g : copy) delete g;
			
		improved = false;
		for (unsigned int i = 0; i < p3.size(); i++){
			p3[i]->evaluate(problem.evalFunc);
			evaluations++;
			if (evaluations >= evalBudget || coco_problem_final_target_hit(problem.PROBLEM)){
				reset();
				for (Genome* g : p0) delete g;
				for (Genome* g : p3) delete g;
				return;
			}

			if (p3[i]->getFitness() < particles[i]->getFitness()){
				particles[i]->setPosition(p3[i]->getX(), p3[i]->getFitness());
			}

			if (p3[i]->getFitness() < bestFitness){
				improved = true;
				bestFitness = p3[i]->getFitness();
			}
		}

		improved ? notImproved = 0 : notImproved++;

		for (Genome* g : p3) delete g;
		for (Genome* g : p0) delete g;

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
