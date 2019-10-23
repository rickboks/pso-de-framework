#pragma once
#include <vector>
#include <fstream>
#include "problem.h"
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include "logger.h"

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
		
		void runSynchronous(Problem problem, int const evalBudget, int const popSize, std::map<int,double> particleUpdateParams);
		void runAsynchronous(Problem problem, int const evalBudget, int const popSize, std::map<int,double> particleUpdateParams);
		public:
		ParticleSwarm(UpdateManagerType const updateManagerType, 
			Topology topologyManager, Synchronicity const synchronous);

		~ParticleSwarm();

		void run(Problem problem, int const evalBudget, int const popSize, std::map<int,double> particleUpdateParams);
		void reset();
		void printParticles();
		void log(std::string filename);
		std::string getIdString() const;
};