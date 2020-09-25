#include <IOHprofiler_problem.h>
#include <IOHprofiler_csv_logger.h>
#include "hybridalgorithm.h"
#include "particle.h"
#include "topologymanager.h"
#include "particleupdatesettings.h"
#include "util.h"
#include "mutationmanager.h"
#include "repairhandler.h"
#include "psode.h"
#include <limits>
#include <iostream>
#include <algorithm> 

PSODE::PSODE(HybridConfig const config)
	: HybridAlgorithm(config){
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

PSODE::~PSODE(){}

void PSODE::run(std::shared_ptr<IOHprofiler_problem<double>> const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger,
    		int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams){
	if (config.synchronicity == "S")
		runSynchronous(problem, logger, evalBudget, popSize, particleUpdateParams);
	else
		runAsynchronous(problem, logger, evalBudget, popSize, particleUpdateParams);
}

void PSODE::runSynchronous(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger, int const evalBudget, 
    		int const popSize, std::map<int,double> const particleUpdateParams){

	std::vector<Particle*> p0;
	std::vector<Particle*> p1;
	std::vector<Particle*> p2;
		
	int const D = problem->IOHprofiler_get_number_of_variables();
	std::vector<double> const lowerBound = problem->IOHprofiler_get_lowerbound();
	std::vector<double> const upperBound = problem->IOHprofiler_get_upperbound();

	ConstraintHandler const*const deCH = deCHs.at(config.deCH)(lowerBound, upperBound);
	ConstraintHandler const*const psoCH = psoCHs.at(config.psoCH)(lowerBound, upperBound);

	ParticleUpdateSettings const settings(config.update, particleUpdateParams, psoCH);

	for (int i = 0; i < popSize; i++){
		particles.push_back(new Particle(D, &settings));
		particles[i]->randomize(lowerBound, upperBound);
	}

	TopologyManager *const topologyManager = topologies.at(config.topology)(particles);
	MutationManager *const mutationManager = mutations.at(config.mutation)(D, deCH);
	CrossoverManager const*const crossoverManager = crossovers.at(config.crossover)(D);
	DEAdaptationManager *const adaptationManager = deAdaptations.at(config.adaptation)();
	SelectionManager const*const selectionManager = selections.at(config.selection)(D, adaptationManager);

	std::vector<double> Fs(popSize);
	std::vector<double> Crs(popSize);	

	while (problem->IOHprofiler_get_evaluations() < evalBudget &&
			!problem->IOHprofiler_hit_optimal()){
		
		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);
		adaptationManager->reset();

		for (int i = 0; i < popSize; i++)
			particles[i]->evaluate(problem,logger);

		p0 = copyPopulation(particles);
		
		for (int i = 0; i < popSize; i++)
			p0[i]->updatePbest();
		for (int i = 0; i < popSize; i++)
			p0[i]->updateGbest();
		for (int i = 0; i < popSize; i++)
			p0[i]->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/double(evalBudget));

		p1 = mutationManager->mutate(particles, Fs);
		p2 = crossoverManager->crossover(particles,p1, Crs);
			
		for (unsigned int i = 0; i < p2.size(); i++){
			p0[i]->evaluate(problem,logger);
			p2[i]->evaluate(problem,logger);
		}

		selectionManager->selection(particles, p0, p2);
		adaptationManager->update();

		for (int i = 0; i < popSize; i++) {
			delete p0[i]; 
			delete p1[i];
			delete p2[i];
		}

		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/double(evalBudget));	
	}

	delete topologyManager;
	delete mutationManager;
	delete crossoverManager;
	delete adaptationManager;
	delete deCH;
	delete psoCH;
	
	for (Particle* const particle : particles)
		delete particle;

	particles.clear();
}

void PSODE::runAsynchronous(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger, int const evalBudget, 
    		int const popSize, std::map<int,double> const particleUpdateParams){

	std::vector<Particle*> p0;
	std::vector<Solution*> p1;
	std::vector<Solution*> p2;

	int const D = problem->IOHprofiler_get_number_of_variables(); /// dimension
	std::vector<double> const lowerBound = problem->IOHprofiler_get_lowerbound();
	std::vector<double> const upperBound = problem->IOHprofiler_get_upperbound();

	ConstraintHandler const*const deCH = deCHs.at(config.deCH)(lowerBound, upperBound);
	ConstraintHandler const*const psoCH = psoCHs.at(config.deCH)(lowerBound, upperBound);
	ParticleUpdateSettings const settings(config.update, particleUpdateParams, psoCH);

	for (int i = 0; i < popSize; i++)
		particles.push_back(new Particle(D, &settings));

	for (Particle* const p : particles)
		p->randomize(lowerBound, upperBound);

	TopologyManager *const topologyManager = topologies.at(config.topology)(particles);
	MutationManager *const mutationManager = mutations.at(config.mutation)(D, deCH);
	CrossoverManager const *const crossoverManager = crossovers.at(config.crossover)(D);
	DEAdaptationManager *const adaptationManager = deAdaptations.at(config.adaptation)();
	SelectionManager const*const selectionManager = selections.at(config.selection)(D, adaptationManager);

	std::vector<double> Fs(popSize);
	std::vector<double> Crs(popSize);

	while (	problem->IOHprofiler_get_evaluations() < evalBudget &&
			!problem->IOHprofiler_hit_optimal()){
		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);
		adaptationManager->reset();
		
		for (int i = 0; i < popSize; i++)			
			particles[i]->evaluate(problem,logger);

		p0 = copyPopulation(particles);
		
		for (int i = 0; i < popSize; i++){
			p0[i]->updatePbest();
			p0[i]->updateGbest();			
			p0[i]->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/double(evalBudget));
		}		

		p1 = mutationManager->mutate(particles, Fs);
		p2 = crossoverManager->crossover(particles,p1, Crs);
			
		for (unsigned int i = 0; i < p2.size(); i++){
			p0[i]->evaluate(problem,logger);
			p2[i]->evaluate(problem,logger);
		}

		selectionManager->selection(particles, p0, p2);
		adaptationManager->update();

		for (int i = 0; i < popSize; i++) {
			delete p0[i]; 
			delete p1[i];
			delete p2[i];
		}

		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/double(evalBudget));	
	}

	delete topologyManager;
	delete mutationManager;
	delete crossoverManager;
	delete adaptationManager;
	delete deCH;
	delete psoCH;
	
	for (Particle* const particle : particles)
		delete particle;

	particles.clear();
}

std::string PSODE::getIdString() const{
	return "H_" + config.update + "_" + config.topology + "_" + config.psoCH + "_" + config.synchronicity + "_" + 
		config.mutation + "_" + config.crossover + "_"  + config.selection + "_" + config.adaptation + "_" + config.deCH; 

}
