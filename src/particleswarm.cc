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
}

void ParticleSwarm::reset(){}

ParticleSwarm::~ParticleSwarm(){}

void ParticleSwarm::run(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger,
    		int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams){

	if (config.synchronicity == "A")
		runAsynchronous(problem, logger, evalBudget, popSize, particleUpdateParams);
	else 
		runSynchronous(problem, logger, evalBudget, popSize, particleUpdateParams);
}

void ParticleSwarm::runAsynchronous(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger,
    		int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams){

	int const D = problem->IOHprofiler_get_number_of_variables(); 
	std::vector<double> lowerBound = problem->IOHprofiler_get_lowerbound(); 
	std::vector<double> upperBound = problem->IOHprofiler_get_upperbound();

	PSOConstraintHandler const* const psoCH = psoCHs.at(config.constraintHandler)(lowerBound, upperBound); 
	ParticleUpdateSettings const settings(config.update, particleUpdateParams, psoCH);

	particles.resize(popSize);
	for (int i = 0; i < popSize; i++){
		particles[i] = new Particle(D, &settings);
		particles[i]->randomize(lowerBound, upperBound);
	}

	TopologyManager* const topologyManager = topologies.at(config.topology)(particles);

	while (	problem->IOHprofiler_get_evaluations() < evalBudget &&
			!problem->IOHprofiler_hit_optimal()){
		
		for (int i = 0; i < popSize; i++){
			particles[i]->evaluate(problem,logger);
			particles[i]->updatePbest();
			particles[i]->updateGbest();
			particles[i]->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/evalBudget);			
		}
		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/evalBudget);	
	}

	delete topologyManager;
	for (Particle* const particle : particles)
		delete particle;
	particles.clear();
	delete psoCH;
}

void ParticleSwarm::runSynchronous(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger,
    		int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams){

	int const D = problem->IOHprofiler_get_number_of_variables(); 
	std::vector<double> const lowerBound = problem->IOHprofiler_get_lowerbound(); 
	std::vector<double> const upperBound = problem->IOHprofiler_get_upperbound(); 

	PSOConstraintHandler* const psoCH = psoCHs.at(config.constraintHandler)(lowerBound, upperBound); 
	ParticleUpdateSettings const settings(config.update, particleUpdateParams, psoCH);

	particles.resize(popSize);
	for (int i = 0; i < popSize; i++){
		particles[i] = new Particle(D, &settings);
		particles[i]->randomize(lowerBound, upperBound);
	}

	TopologyManager* const topologyManager = topologies.at(config.topology)(particles);

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
	
		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/evalBudget);	
	}

	delete topologyManager;
	for (Particle* const particle : particles)
		delete particle;
	particles.clear();
	delete psoCH;
}

std::string ParticleSwarm::getIdString() const {
	return "P_" + config.update + "_" + config.topology + "_" + config.constraintHandler + "_" + config.synchronicity;
}
