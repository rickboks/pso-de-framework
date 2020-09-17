#pragma once
#include <vector>
#include <fstream>
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include <memory>

struct Problem;
class Particle;
struct ParticleUpdateSettings;
template <typename T> 
class IOHprofiler_problem;
class IOHprofiler_csv_logger;

struct PSOConfig {
	PSOConfig(std::string const update, std::string const topology, std::string constraintHandler, bool const synchronous) 
	: update(update), topology(topology), constraintHandler(constraintHandler), synchronous(synchronous){}
	std::string const update, topology, constraintHandler;
	bool const synchronous;
};

class ParticleSwarm {
	private:
		PSOConfig const config;
		std::vector<Particle*> particles;
		TopologyManager* topologyManager;
		ConstraintHandler* psoCH;

		std::shared_ptr<IOHprofiler_problem<double> > problem;
    	std::shared_ptr<IOHprofiler_csv_logger> logger;
		bool logging;

		void runSynchronous(int const evalBudget, int const popSize, 
    		std::map<int,double> particleUpdateParams);
		void runAsynchronous(int const evalBudget, int const popSize, 
    		std::map<int,double> particleUpdateParams);
	public:
		ParticleSwarm(PSOConfig config);
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
