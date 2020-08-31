#pragma once
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "selectionmanager.h"
#include <memory>

class HybridAlgorithm {
	protected:
		UpdateManagerType const updateManagerType;
		Topology topologyManagerType;
		std::vector<Particle*> particles;		
		Synchronicity const synchronicity;
		MutationType const mutationType;
		CrossoverType const crossoverType;
		SelectionType const selectionType;
		DEAdaptationType const adaptationType;
		TopologyManager* topologyManager;
		MutationManager<Particle>* mutationManager;
		CrossoverManager<Particle>* crossoverManager;
		DEAdaptationManager* adaptationManager;
		SelectionManager* selectionManager;
		std::shared_ptr<IOHprofiler_problem<double> > problem; 
		std::shared_ptr<IOHprofiler_csv_logger> logger;
		public:
		HybridAlgorithm(UpdateManagerType const updateManagerType, 
			Topology topologyManager, Synchronicity const synchronous, 
			MutationType const mutationType, CrossoverType const crossoverType, 
			SelectionType const selection, DEAdaptationType adaptionType);

		virtual ~HybridAlgorithm() = 0;

		virtual void run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger, int const evalBudget, 
    		int const popSize, std::map<int,double> particleUpdateParams) = 0;

		virtual void reset() = 0;
		virtual std::string getIdString() const = 0;
};
