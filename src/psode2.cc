#include <IOHprofiler_problem.h>
#include <IOHprofiler_csv_logger.h>
#include "util.h"
#include "constrainthandler.h"
#include "hybridalgorithm.h"
#include "repairhandler.h"
#include "particle.h"
#include "topologymanager.h"
#include "particleupdatesettings.h"
#include "util.h"
#include "mutationmanager.h"
#include "psode2.h"
#include "deadaptationmanager.h"
#include <limits>
#include <iostream>
#include <algorithm> 

PSODE2::PSODE2(HybridConfig const config)
		: HybridAlgorithm(config){}

PSODE2::~PSODE2(){}

void PSODE2::run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger,
    		int const evalBudget, int const popSize, std::map<int,double> particleUpdateParams){
	runAsynchronous(problem, logger, evalBudget, popSize, particleUpdateParams);
}

void PSODE2::runAsynchronous(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger, int const evalBudget, int const popSize, 
			std::map<int,double> particleUpdateParams){

	int const D = problem->IOHprofiler_get_number_of_variables(); /// dimension

	std::vector<double> const lowerBound = problem->IOHprofiler_get_lowerbound();
	std::vector<double> const upperBound = problem->IOHprofiler_get_upperbound();

	DEConstraintHandler *const deCH = deCHs.at(config.deCH)(lowerBound,upperBound);
	PSOConstraintHandler *const psoCH = psoCHs.at(config.psoCH)(lowerBound,upperBound);
	ParticleUpdateSettings const settings(config.update, particleUpdateParams, psoCH);

	int const split = popSize / 2;
	for (int i = 0; i < split; i++) psoPop.push_back(new Particle(D, &settings));
	for (int i = split; i < popSize; i++) dePop.push_back(new Solution(D));
	particles.insert(particles.end(), psoPop.begin(), psoPop.end());
	particles.insert(particles.end(), dePop.begin(), dePop.end());

	for (Solution* const p : particles){
		p->randomize(lowerBound, upperBound);
		p->evaluate(problem, logger);
	}

	TopologyManager* const topologyManager = topologies.at(config.topology)(psoPop);
	MutationManager* const mutationManager = mutations.at(config.mutation)(D, deCH);
	CrossoverManager const*const crossoverManager = crossovers.at(config.crossover)(D);
	DEAdaptationManager *const adaptationManager = deAdaptations.at(config.adaptation)(popSize);

	std::vector<double> Fs(dePop.size());
	std::vector<double> Crs(dePop.size());

	int iterations = 0;
	while (problem->IOHprofiler_get_evaluations() < evalBudget &&
			!problem->IOHprofiler_hit_optimal()){

		// Get new DE parameters from the adaptation manager (JADE or constant)
		adaptationManager->nextF(Fs);
		adaptationManager->nextCr(Crs);

		// Measure fitness of the PSO population
		// Update pbest and gbest of PSO population
		// Update velocity and position of PSO population
		for (Particle* const p : psoPop){
			p->updatePbest();
			p->updateGbest();
			p->updateVelocityAndPosition(double(problem->IOHprofiler_get_evaluations())/double(evalBudget));			
			p->evaluate(problem,logger);
		}

		// Perform mutation 
		std::vector<Solution*> const donors = mutationManager->mutate(dePop, Fs);
		// Perform crossover
		std::vector<Solution*> const trials = crossoverManager->crossover(dePop, donors, Crs);

		for (Solution* d : donors) 
			delete d;

		std::vector<double> parentF(popSize), trialF(popSize);
		for (unsigned int i = 0; i < dePop.size(); i++){
			//Evaluate the parent vector
			parentF[i] = dePop[i]->evaluate(problem,logger);

			//Evaluate the trial vector
			trialF[i] = trials[i]->evaluate(problem,logger);

			// Perform selection
			if ( trialF[i] < parentF[i] ){
				dePop[i]->setX(trials[i]->getX(), trials[i]->getFitness());
			}
		}

		for (Solution* p : trials)
			delete p;

		if (iterations % 10 == 0)
			share();

		adaptationManager->update(parentF, trialF);
		iterations++;	
		topologyManager->update(double(problem->IOHprofiler_get_evaluations())/evalBudget);	
	}

	delete topologyManager;
	delete mutationManager;
	delete crossoverManager;
	delete adaptationManager;
	delete deCH;
	delete psoCH;

	for (Solution* particle : particles)
		delete particle;

	particles.clear();
	dePop.clear();
	psoPop.clear();
}

std::string PSODE2::getIdString() const{
	return "H2_" + config.update + "_" + config.topology + "_" + config.psoCH /*+ "_" + config.synchronicity*/
		+ "_" + config.mutation + "_" + config.crossover + "_" + config.adaptation + "_" + config.deCH;
}

void PSODE2::share(){
	Solution* const best_de = getPBest(dePop);
	Particle* const best_pso = getPBest(psoPop);

	std::vector<double> const x = best_pso->getX();
	double const y = best_pso->getFitness();

	best_pso->setXandUpdateV(best_de->getX(), best_de->getFitness());
	best_de->setX(x, y);
}
