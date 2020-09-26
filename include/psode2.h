#pragma once
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
//#include "selectionmanager.h"
#include "hybridalgorithm.h"
#include <memory>

class PSODE2 : public HybridAlgorithm {
	private:
		std::vector<Solution*> particles;		
		std::vector<Particle*> psoPop;
		std::vector<Solution*> dePop;

		void runAsynchronous(std::shared_ptr<IOHprofiler_problem<double>> const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger, int const evalBudget, 
    		int const popSize, std::map<int,double> const particleUpdateParams);
		void share();

	public:
		PSODE2(HybridConfig const config);
		~PSODE2();

		void run(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger, int const evalBudget, 
    		int const popSize, std::map<int,double> const particleUpdateParams);

		void reset();
		std::string getIdString() const;
};
