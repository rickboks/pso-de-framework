#pragma once
#include <vector>
#include <map>
#include "particleupdatemanager.h"
#include "topologymanager.h"
#include "particleswarm.h"

class ParticleSwarmSuite {
	private:		
		std::vector<std::string> updateManagers;
		std::vector<std::string> topologyManagers;
		std::vector<std::string> constraintHandlers;
		std::vector<std::string> synchronicities;
		std::vector<PSOConfig> configurations;
	public:
		ParticleSwarmSuite();
		void setUpdateManagers(std::vector<std::string> updateManagers);
		void setTopologyManagers(std::vector<std::string> topologyManagers);
		void setSynchronicities(std::vector<std::string> synchronicities);
		void setConstraintHandlers(std::vector<std::string> chs);
		void generateConfigurations();
		ParticleSwarm getParticleSwarm(int const i);	
		int size() const;
};
