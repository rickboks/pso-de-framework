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
	PSOConfig(std::string const update, std::string const topology, std::string constraintHandler, std::string const synchronicity) 
	: update(update), topology(topology), constraintHandler(constraintHandler), synchronicity(synchronicity){}

	std::string const update, topology, constraintHandler, synchronicity;
};

class ParticleSwarm {
	private:
		PSOConfig const config;

		void runSynchronous(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger,
    		int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams);

		void runAsynchronous(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger,
    		int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams);
	public:
		ParticleSwarm(PSOConfig const config);
		~ParticleSwarm();

		void run(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
    		std::shared_ptr<IOHprofiler_csv_logger> const logger,
    		int const evalBudget, int const popSize, std::map<int,double> const particleUpdateParams);

		void reset();
		std::string getIdString() const;
};
