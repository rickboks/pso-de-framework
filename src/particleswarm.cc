#include "particleswarm.h"
#include "particle.h"
#include "searchspace.h"
#include "topologymanager.h"
#include "particleupdatesettings.h"
#include "coco.h"
#include <limits>
#include <iostream>
#include <problem.h>

constexpr double DOUBLE_MIN = std::numeric_limits<double>::min();
constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();

ParticleSwarm::ParticleSwarm(UpdateManagerType const updateManagerType,
	Topology topologyManagerType, Synchronicity const synchronicity):

	updateManagerType(updateManagerType),
	topologyManagerType(topologyManagerType), topologyManager(NULL), synchronicity(synchronicity){
}

void ParticleSwarm::reset(){
	delete topologyManager;
	topologyManager = NULL;
	
	for (Particle* const particle : particles)
		delete particle;

	particles.clear();
}

ParticleSwarm::~ParticleSwarm(){
	if (topologyManager != NULL)
		delete topologyManager;

	if (!particles.empty())
	for (Particle* const particle : particles)
		delete particle;
}

void ParticleSwarm::run(Problem const problem, int const evalBudget, int popSize, std::map<int,double> particleUpdateParams){
	if (synchronicity == SYNCHRONOUS)
		runSynchronous(problem, evalBudget, popSize, particleUpdateParams);
	else if (synchronicity == ASYNCHRONOUS)
		runAsynchronous(problem, evalBudget, popSize, particleUpdateParams);
}

void ParticleSwarm::runAsynchronous(Problem const problem, int const evalBudget, int popSize, std::map<int,double> particleUpdateParams){
	this->topologyManager = TopologyManagerFactory::createTopologyManager(topologyManagerType, particles);
	popSize = this->topologyManager->getClosestValidPopulationSize(popSize);

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
			std::vector<double> position = particles[i]->getPosition();			

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

			particles[i]->updatePbest();
			particles[i]->updateGbest();
			particles[i]->updateVelocityAndPosition(double(evaluations)/double(evalBudget));			
		}

		if (improved)
			notImproved = 0;
		else 
			notImproved++;	

		iterations++;	
		topologyManager->update(double(evaluations)/double(evalBudget));	
	}

	reset();
}

void ParticleSwarm::runSynchronous(Problem const problem, int const evalBudget, int popSize, std::map<int,double> particleUpdateParams){
		this->topologyManager = TopologyManagerFactory::createTopologyManager(topologyManagerType, particles);
	popSize = this->topologyManager->getClosestValidPopulationSize(popSize);

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
			std::vector<double> position = particles[i]->getPosition();			

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

		if (improved)
			notImproved = 0;
		else 
			notImproved++;
		
		for (int i = 0; i < popSize; i++)
			particles[i]->updatePbest();

		for (int i = 0; i < popSize; i++)
			particles[i]->updateGbest();

		for (int i = 0; i < popSize; i++)
			particles[i]->updateVelocityAndPosition(double(evaluations)/double(evalBudget));
	
		iterations++;	
		topologyManager->update(double(evaluations)/double(evalBudget));	
	}

	reset();
}

std::string ParticleSwarm::getIdString() const {
	std::string id = "PS_";

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

	if(synchronicity == SYNCHRONOUS)
		id+= "S";
	else
		id+= "A";

	return id;
}