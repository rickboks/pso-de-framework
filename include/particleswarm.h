#pragma once
#include <vector>
#include <fstream>
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include <memory>


enum Synchronicity {
	SYNCHRONOUS,
	ASYNCHRONOUS,
	SYNC_END
};

struct Problem;
class Particle;
struct ParticleUpdateSettings;
template <typename T> 
class IOHprofiler_problem;
class IOHprofiler_csv_logger;

class ParticleSwarm {
	private:
		UpdateManagerType const updateManagerType;
		Topology topologyManagerType;
		std::vector<Particle*> particles;
		TopologyManager* topologyManager;
		Synchronicity const synchronicity;
		std::shared_ptr<IOHprofiler_problem<double> > problem;
    	std::shared_ptr<IOHprofiler_csv_logger> logger;
		bool logging;

		void runSynchronous(int const evalBudget, int const popSize, 
    		std::map<int,double> particleUpdateParams);
		void runAsynchronous(int const evalBudget, int const popSize, 
    		std::map<int,double> particleUpdateParams);
	public:
		ParticleSwarm(UpdateManagerType const updateManagerType, 
			Topology topologyManager, Synchronicity const synchronous);

		~ParticleSwarm();

		void run(std::shared_ptr<IOHprofiler_problem<double> > problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> logger,
    		int const evalBudget, int const popSize, std::map<int,double> particleUpdateParams);

		void reset();
		std::string getIdString() const;
		void enableLogging();
		void logStart();
		void logEnd();
		void logPositions();
};
