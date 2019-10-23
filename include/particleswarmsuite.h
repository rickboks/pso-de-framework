#pragma once

#include <vector>
#include <map>
#include "particleupdatemanager.h"
#include "topologymanager.h"
#include "particleswarm.h"
typedef std::tuple<UpdateManagerType, Topology, Synchronicity> pso_configuration;

class ParticleSwarmSuite {
	private:		
		std::vector<UpdateManagerType> updateManagers;
		std::vector<Topology> topologyManagers;
		std::vector<Synchronicity> synchronicities;

		std::vector<pso_configuration> configurations;
		public:
		ParticleSwarmSuite();
		void setUpdateManagers(std::vector<UpdateManagerType> updateManagers);
		void setTopologyManagers(std::vector<Topology> topologyManagers);
		void setSynchronicities(std::vector<Synchronicity> synchronicities);
		ParticleSwarm getParticleSwarm(int const i);	
		int size() const;
};