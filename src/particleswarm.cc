#include <IOHprofiler_problem.h>
#include <IOHprofiler_csv_logger.h>
#include "particleswarm.h"
#include "particle.h"
#include "topologymanager.h"
#include "particleupdatesettings.h"
#include "repairhandler.h"
#include <limits>
#include <iostream>
#include <fstream>

ParticleSwarm::ParticleSwarm(PSOConfig const config) : config(config), topologyManager(NULL){
	logging = false;
}

void ParticleSwarm::reset(){
	delete topologyManager;
	topologyManager = NULL;
	
	for (Particle* const particle : particles)
		delete particle;

	particles.clear();

	delete psoCH;
}

ParticleSwarm::~ParticleSwarm(){
	if (topologyManager != NULL)
		delete topologyManager;

	if (!particles.empty())
		for (Particle* const particle : particles)
			delete particle;
}

void ParticleSwarm::run(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger,
    		int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams){

	this->problem=problem;
	this->logger=logger;

	if (config.synchronicity == "A")
		runSynchronous(evalBudget, popSize, particleUpdateParams);
	else 
		runAsynchronous(evalBudget, popSize, particleUpdateParams);
}

void ParticleSwarm::runAsynchronous(int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams){
	int const D = problem->IOHprofiler_get_number_of_variables(); 
	std::vector<double> lowerBound = problem->IOHprofiler_get_lowerbound(); 
	std::vector<double> upperBound = problem->IOHprofiler_get_upperbound();

	psoCH = psoCHs.at(config.constraintHandler)(lowerBound, upperBound); 
	ParticleUpdateSettings settings(config.update, particleUpdateParams, psoCH);

	for (int i = 0; i < popSize; i++){
		Particle* p = new Particle(D, &settings);
		p->randomize(lowerBound, upperBound);
		particles.push_back(p);
	}

	logStart();
	logPositions();

	topologyManager = topologies.at(config.topology)(particles);
	
	int iterations = 0;
	double bestFitness = std::numeric_limits<double>::max();
	int notImproved = 0;
	bool improved;

	while (	problem->IOHprofiler_get_evaluations() < evalBudget &&
			!problem->IOHprofiler_hit_optimal()){
		
		improved = false;

		for (int i = 0; i < popSize; i++){
			double const y = particles[i]->evaluate(problem,logger);

			if (y < bestFitness){
				improved = true;
				bestFitness = y;
			}

			particles[i]->updatePbest();
			particles[i]->updateGbest();
			particles[i]->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/evalBudget);			
		}
		logPositions();

		improved ? notImproved=0 : notImproved++;

		iterations++;	
		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/evalBudget);	
	}

	reset();
	logEnd();
}

void ParticleSwarm::runSynchronous(int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams){
	int const D = problem->IOHprofiler_get_number_of_variables(); 
	std::vector<double> lowerBound = problem->IOHprofiler_get_lowerbound(); 
	std::vector<double> upperBound = problem->IOHprofiler_get_upperbound(); 

	ConstraintHandler* psoCH = psoCHs.at(config.constraintHandler)(lowerBound, upperBound); 
	ParticleUpdateSettings settings(config.update, particleUpdateParams, psoCH);

	for (int i = 0; i < popSize; i++){
		Particle* p = new Particle(D, &settings);
		p->randomize(lowerBound, upperBound);
		particles.push_back(p);
	}

	topologyManager = topologies.at(config.topology)(particles);

	logStart();
	logPositions();

	int iterations = 0;
	double bestFitness = std::numeric_limits<double>::max();
	int notImproved = 0;
	bool improved;

	while (	problem->IOHprofiler_get_evaluations() < evalBudget &&
			!problem->IOHprofiler_hit_optimal()){
		
		improved = false;

		for (int i = 0; i < popSize; i++){		
			double const y = particles[i]->evaluate(problem,logger);

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
			particles[i]->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/evalBudget);
	
		logPositions();
	
		iterations++;	
		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/evalBudget);	
	}

	reset();
	logEnd();
}

std::string ParticleSwarm::getIdString() const {
	return "P_" + config.update + "_" + config.topology + "_" + config.constraintHandler + "_" + config.synchronicity;
}

void ParticleSwarm::enableLogging(){
	logging=true;
}

void ParticleSwarm::logStart(){
	if (logging)
		std::cout << "LOGGER: START" << std::endl;
}

void ParticleSwarm::logEnd(){
	if (logging)
		std::cout << "LOGGER: END" << std::endl;
}

void ParticleSwarm::logPositions(){
	if (logging){
		std::cout << "LOGGER: " << std::endl;
		for (unsigned int i = 0; i < particles.size(); i++){
			std::cout << "LOGGER: ";
			std::vector<double> position = particles[i]->getX();
			for (unsigned int j = 0; j < position.size()-1 ; j++){
				std::cout << position[j] << " ";
			}
			std::cout << position[position.size()-1] << std::endl;
		}
	}
}
