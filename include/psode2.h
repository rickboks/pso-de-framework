#pragma once
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "selectionmanager.h"
#include "hybridalgorithm.h"
#include <memory>

class PSODE2 : public HybridAlgorithm {
	private:
		std::vector<Particle*> particles;		
		//void runSynchronous(int const evalBudget, int const popSize, 
			//std::map<int,double> particleUpdateParams);

		void runAsynchronous(int const evalBudget, int const popSize, 
			std::map<int,double> particleUpdateParams);

		std::vector<Particle*> copyPopulation(std::vector<Particle*>const& particles);

		public:
		PSODE2(UpdateManagerType const updateManagerType, 
			Topology topologyManager, Synchronicity const synchronous, 
			MutationType const mutationType, CrossoverType const crossoverType, 
			SelectionType const selection, DEAdaptationType adaptionType);

		PSODE2(hybrid_config config);

		~PSODE2();

		void run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger, int const evalBudget, 
    		int const popSize, std::map<int,double> particleUpdateParams);

		void reset();
		std::string getIdString() const;
};
