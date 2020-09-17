#pragma once
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "selectionmanager.h"
#include "constrainthandler.h"
#include <memory>

struct HybridConfig {
	std::string const update, topology, psoCH;
	bool const synchronous;
	std::string const  mutation, crossover, selection, adaptation, deCH;

	HybridConfig(std::string const update, std::string const topology, std::string const psoCH, bool const synchronous, 
			std::string const mutation, std::string const crossover, std::string const selection, 
			std::string const adaptation, std::string const deCH):
				update(update), topology(topology), psoCH(psoCH), synchronous(synchronous),
				mutation(mutation), crossover(crossover), selection(selection),
				adaptation(adaptation), deCH(deCH){}
};

class HybridAlgorithm {
	protected:
		HybridConfig const config;
		TopologyManager* topologyManager;
		MutationManager* mutationManager;
		CrossoverManager* crossoverManager;
		DEAdaptationManager* adaptationManager;
		SelectionManager* selectionManager;
		ConstraintHandler* deCH;
		ConstraintHandler* psoCH;
		std::shared_ptr<IOHprofiler_problem<double> > problem; 
		std::shared_ptr<IOHprofiler_csv_logger> logger;

	public:
		HybridAlgorithm(HybridConfig config);
		virtual ~HybridAlgorithm() = 0;

		virtual void run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
		std::shared_ptr<IOHprofiler_csv_logger> logger, int const evalBudget, 
		int const popSize, std::map<int,double> particleUpdateParams) = 0;

		virtual void reset() = 0;
		virtual std::string getIdString() const = 0;
};
