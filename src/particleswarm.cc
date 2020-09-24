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

ParticleSwarm::ParticleSwarm(PSOConfig const config) : config(config){
	logging = false;
}

void ParticleSwarm::reset(){}

ParticleSwarm::~ParticleSwarm(){}

void ParticleSwarm::run(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger,
    		int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams){

	if (config.synchronicity == "A")
		runSynchronous(problem, logger, evalBudget, popSize, particleUpdateParams);
	else 
		runAsynchronous(problem, logger, evalBudget, popSize, particleUpdateParams);
}

void ParticleSwarm::runAsynchronous(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger,
    		int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams){

	int const D = problem->IOHprofiler_get_number_of_variables(); 
	std::vector<double> lowerBound = problem->IOHprofiler_get_lowerbound(); 
	std::vector<double> upperBound = problem->IOHprofiler_get_upperbound();

	ConstraintHandler const* const psoCH = psoCHs.at(config.constraintHandler)(lowerBound, upperBound); 
	ParticleUpdateSettings const settings(config.update, particleUpdateParams, psoCH);

	for (int i = 0; i < popSize; i++){
		Particle* p = new Particle(D, &settings);
		p->randomize(lowerBound, upperBound);
		particles.push_back(p);
	}

	logStart();
	logPositions();

	TopologyManager* const topologyManager = topologies.at(config.topology)(particles);
	
	while (	problem->IOHprofiler_get_evaluations() < evalBudget &&
			!problem->IOHprofiler_hit_optimal()){
		
		for (int i = 0; i < popSize; i++){
			particles[i]->evaluate(problem,logger);
			particles[i]->updatePbest();
			particles[i]->updateGbest();
			particles[i]->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/evalBudget);			
		}

		logPositions();
		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/evalBudget);	
	}

	delete topologyManager;
	for (Particle* const particle : particles)
		delete particle;
	particles.clear();
	delete psoCH;

	logEnd();
}

void ParticleSwarm::runSynchronous(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger,
    		int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams){

	int const D = problem->IOHprofiler_get_number_of_variables(); 
	std::vector<double> const lowerBound = problem->IOHprofiler_get_lowerbound(); 
	std::vector<double> const upperBound = problem->IOHprofiler_get_upperbound(); 

	ConstraintHandler* const psoCH = psoCHs.at(config.constraintHandler)(lowerBound, upperBound); 
	ParticleUpdateSettings const settings(config.update, particleUpdateParams, psoCH);

	for (int i = 0; i < popSize; i++){
		Particle* p = new Particle(D, &settings);
		p->randomize(lowerBound, upperBound);
		particles.push_back(p);
	}

	TopologyManager* const topologyManager = topologies.at(config.topology)(particles);

	logStart();
	logPositions();

	while (	problem->IOHprofiler_get_evaluations() < evalBudget &&
			!problem->IOHprofiler_hit_optimal()){
		
		for (int i = 0; i < popSize; i++){
			particles[i]->evaluate(problem,logger);
			particles[i]->updatePbest();
		}

		for (int i = 0; i < popSize; i++){
			particles[i]->updateGbest();
			particles[i]->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/evalBudget);
		}
	
		logPositions();
		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/evalBudget);	
	}

	delete topologyManager;
	for (Particle* const particle : particles)
		delete particle;
	particles.clear();
	delete psoCH;

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
