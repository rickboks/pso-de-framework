#include "particleswarm.h"
#include "particle.h"
#include "topologymanager.h"
#include "particleupdatesettings.h"
#include "instancenamer.h"
#include "coco.h"
#include <limits>
#include <iostream>
#include <problem.h>
#include <fstream>

constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();

ParticleSwarm::ParticleSwarm(UpdateManagerType const updateManagerType,
	Topology topologyManagerType, Synchronicity const synchronicity):

	updateManagerType(updateManagerType),
	topologyManagerType(topologyManagerType), topologyManager(NULL), synchronicity(synchronicity), 
	outfile(new std::ofstream()){
}

void ParticleSwarm::reset(){
	delete topologyManager;
	topologyManager = NULL;
	
	for (Particle* const particle : particles)
		delete particle;

	particles.clear();
}

ParticleSwarm::~ParticleSwarm(){

	if (outfile->is_open())
		outfile->close();
	delete outfile;
	outfile = NULL;
	
	if (topologyManager != NULL)
		delete topologyManager;

	if (!particles.empty())
		for (Particle* const particle : particles)
			delete particle;

}

void ParticleSwarm::run(Problem problem, int const evalBudget, int popSize, std::map<int,double> particleUpdateParams){
	if (synchronicity == SYNCHRONOUS)
		runSynchronous(problem, evalBudget, popSize, particleUpdateParams);
	else if (synchronicity == ASYNCHRONOUS)
		runAsynchronous(problem, evalBudget, popSize, particleUpdateParams);
}

void ParticleSwarm::runAsynchronous(Problem problem, int const evalBudget, 
	int popSize, std::map<int,double> particleUpdateParams){
	this->topologyManager = TopologyManagerFactory::createTopologyManager(topologyManagerType, particles);
	popSize = this->topologyManager->getClosestValidPopulationSize(popSize);

	int const dimension = coco_problem_get_dimension(problem.PROBLEM);
	double const* smallest = coco_problem_get_smallest_values_of_interest(problem.PROBLEM);
	double const* largest = coco_problem_get_largest_values_of_interest(problem.PROBLEM);

	ParticleUpdateSettings settings(updateManagerType, particleUpdateParams, 
				std::vector<double>(smallest, smallest + dimension), 
				std::vector<double>(largest, largest + dimension));

	for (int i = 0; i < popSize; i++){
		Particle* p = new Particle(coco_problem_get_dimension(problem.PROBLEM), settings);
		p->randomize(settings.xMax, settings.xMin);
		particles.push_back(p);
	}
	for (Particle* const p : particles)
		p->updateGbest();

	topologyManager->initialize();
	
	int iterations = 0;
	double bestFitness = DOUBLE_MAX;
	int notImproved = 0;
	bool improved;
	int evaluations = 0;

	while (	
			notImproved < 100 && 
			evaluations <= evalBudget &&
			!coco_problem_final_target_hit(problem.PROBLEM)){

		printParticles();
		
		improved = false;

		for (int i = 0; i < popSize; i++){
			std::vector<double> position = particles[i]->getPosition();			

			double y = particles[i]->evaluate(problem.evalFunc);
			evaluations++;

			// if (evaluations >= evalBudget || coco_problem_final_target_hit(problem.PROBLEM)){
			// 	reset();
			// 	return;
			// }

			if (y < bestFitness){
				improved = true;
				bestFitness = y;
			}

			particles[i]->updatePbest();
			particles[i]->updateGbest();
			particles[i]->updateVelocityAndPosition(double(evaluations)/double(evalBudget));			
		}

		improved ? notImproved=0 : notImproved++;

		iterations++;	
		topologyManager->update(double(evaluations)/double(evalBudget));	
	}

	reset();
}

void ParticleSwarm::runSynchronous(Problem problem, int const evalBudget, int popSize, 
	std::map<int,double> particleUpdateParams){

	this->topologyManager = TopologyManagerFactory::createTopologyManager(topologyManagerType, particles);
	popSize = this->topologyManager->getClosestValidPopulationSize(popSize);

	int const dimension = coco_problem_get_dimension(problem.PROBLEM);
	double const *smallest = coco_problem_get_smallest_values_of_interest(problem.PROBLEM);
	double const *largest = coco_problem_get_largest_values_of_interest(problem.PROBLEM);

	ParticleUpdateSettings settings(updateManagerType, particleUpdateParams, 
				std::vector<double>(smallest, smallest + dimension), 
				std::vector<double>(largest, largest + dimension));

	for (int i = 0; i < popSize; i++){
		Particle* p = new Particle(coco_problem_get_dimension(problem.PROBLEM), settings);
		p->randomize(settings.xMax, settings.xMin);
		particles.push_back(p);
	}
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

			// if (evaluations >= evalBudget || coco_problem_final_target_hit(problem.PROBLEM)){
			// 	reset();
			// 	return;
			// }

			if (y < bestFitness){
				improved = true;
				bestFitness = y;
			}
	
		}		

		improved ? notImproved=0 : notImproved++;
		
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
	return InstanceNamer::getName(updateManagerType, topologyManagerType, synchronicity);
}

void ParticleSwarm::log(std::string file){
	outfile->open(file);
}

void ParticleSwarm::printParticles(){
	if (outfile != NULL && outfile->is_open()){
		for(unsigned int i = 0; i < particles.size(); i++){
			*outfile << particles[i]->positionString() << std::endl;
		}

		*outfile << std::endl;
	}
}