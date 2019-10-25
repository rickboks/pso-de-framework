#pragma once
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "selectionmanager.h"
#include <memory>

class HybridAlgorithm {
	private:
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
		void runSynchronous(int const evalBudget, int const popSize, 
			std::map<int,double> particleUpdateParams,
			double const F, double const Cr);

		void runAsynchronous(int const evalBudget, int const popSize, 
			std::map<int,double> particleUpdateParams,
			double const F, double const Cr);

		public:
		HybridAlgorithm(UpdateManagerType const updateManagerType, 
			Topology topologyManager, Synchronicity const synchronous, 
			MutationType const mutationType, CrossoverType const crossoverType, 
			SelectionType const selection, DEAdaptationType adaptionType);

		~HybridAlgorithm();

		void run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger, int const evalBudget, 
    		int const popSize, std::map<int,double> particleUpdateParams, 
			double const F, double const Cr);

		void reset();
		std::string getIdString() const;

		std::vector<Particle*> copyPopulation(std::vector<Particle*>& particles);

		template <class T>
		std::vector<Solution*> toSolutions(std::vector<T*> vec);
		std::vector<Particle*> toParticles(std::vector<Solution*> sol);
};