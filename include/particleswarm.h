#pragma once
#include <vector>
#include <fstream>
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include <memory>
#include "iohsrc/Template/IOHprofiler_problem.hpp"
#include "iohsrc/Template/Loggers/IOHprofiler_csv_logger.h"

enum Synchronicity {
	SYNCHRONOUS,
	ASYNCHRONOUS,
	SYNC_END
};

struct Problem;
class Particle;
class ParticleUpdateSettings;
class ParticleSwarm {
	private:
		UpdateManagerType const updateManagerType;
		Topology topologyManagerType;
		std::vector<Particle*> particles;
		TopologyManager* topologyManager;
		Synchronicity const synchronicity;
		std::ofstream* outfile;
		
		void runSynchronous(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger, int const evalBudget, int const popSize, 
    		std::map<int,double> particleUpdateParams);
		void runAsynchronous(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger, int const evalBudget, int const popSize, 
    		std::map<int,double> particleUpdateParams);
		public:
		ParticleSwarm(UpdateManagerType const updateManagerType, 
			Topology topologyManager, Synchronicity const synchronous);

		~ParticleSwarm();

		void run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger,
    		int const evalBudget, int const popSize, std::map<int,double> particleUpdateParams);

		void reset();
		void printParticles();
		void log(std::string filename);
		std::string getIdString() const;
};