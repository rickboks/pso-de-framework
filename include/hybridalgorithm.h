#pragma once
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "selectionmanager.h"
#include "constrainthandler.h"
#include <memory>

//typedef std::tuple<UpdateManagerType, Topology, Synchronicity, MutationType, 
	//CrossoverType, SelectionType, DEAdaptationType> hybrid_config;
struct hybrid_config {
	UpdateManagerType const update;
	Topology const topology;
	Synchronicity const synchronicity;
	MutationType const mutation;
	CrossoverType const crossover;
	SelectionType const selection;
	DEAdaptationType const adaptation;
	std::string const psoCH, deCH;

	hybrid_config(UpdateManagerType update, Topology topology,
		Synchronicity synchronicity, MutationType mutation,
		CrossoverType crossover, SelectionType selection, DEAdaptationType adaptation,
		std::string psoCH, std::string deCH):
		update(update), topology(topology), synchronicity(synchronicity),
		mutation(mutation), crossover(crossover), selection(selection),
		adaptation(adaptation), psoCH(psoCH), deCH(deCH){
	}
};

class HybridAlgorithm {
	protected:
		hybrid_config const config;
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
		HybridAlgorithm(UpdateManagerType const updateManagerType, 
			Topology topologyManager, Synchronicity const synchronous, 
			MutationType const mutationType, CrossoverType const crossoverType, 
			SelectionType const selection, DEAdaptationType const adaptionType, std::string const psoCH, std::string const deCH);

		HybridAlgorithm(hybrid_config config);

		virtual ~HybridAlgorithm() = 0;

		virtual void run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
		std::shared_ptr<IOHprofiler_csv_logger> logger, int const evalBudget, 
		int const popSize, std::map<int,double> particleUpdateParams) = 0;

		virtual void reset() = 0;
		virtual std::string getIdString() const = 0;
};
