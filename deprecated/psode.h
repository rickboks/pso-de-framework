#pragma once
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "selectionmanager.h"
#include "hybridalgorithm.h"
#include <memory>

class PSODE : public HybridAlgorithm {
	private:
		std::vector<Particle*> particles;		

		void runSynchronous(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger, int const evalBudget, 
    		int const popSize, std::map<int,double> const particleUpdateParams);

		void runAsynchronous(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger, int const evalBudget, 
    		int const popSize, std::map<int,double> const particleUpdateParams);

		std::vector<Particle*> copyPopulation(std::vector<Particle*>const& particles);

		public:
		PSODE(HybridConfig const config);
		~PSODE();

		void run(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger, int const evalBudget, 
    		int const popSize, std::map<int,double> const particleUpdateParams);

		std::string getIdString() const;
};
