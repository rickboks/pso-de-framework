#pragma once
#include <vector>
#include "searchspace.h"
#include "problem.h"
#include "particleupdatesettings.h"
#include "topologymanager.h"

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
		void runSynchronous(Problem const problem, int const evalBudget, int const popSize, std::map<int,double> particleUpdateParams);
		void runAsynchronous(Problem const problem, int const evalBudget, int const popSize, std::map<int,double> particleUpdateParams);

		public:
		ParticleSwarm(UpdateManagerType const updateManagerType, 
			Topology topologyManager, Synchronicity const synchronous);

		~ParticleSwarm();

		void run(Problem const problem, int const evalBudget, int const popSize, std::map<int,double> particleUpdateParams);
		void reset();
		
		std::string getIdString() const;
};