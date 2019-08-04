#pragma once
#include <vector>
#include <map>
#include "particleswarm.h"
#include "particleupdatemanager.h"
#include "topologymanager.h"

class ParticleSwarmSuite {
	private:
		std::vector<UpdateManagerType> updateManagers;
		std::vector<Topology> topologyManagers;
		std::vector<Synchronicity> synchronicities;
	public:
		ParticleSwarmSuite();

		class ParticleSwarmSuiteIterator {
			private:
				std::vector<UpdateManagerType> const updateManagers;
				std::vector<Topology> const topologyManagers;
				std::vector<Synchronicity> const synchronicities;
				int updateManager;
				int topology;
				int synchronicity;
			public:
				ParticleSwarmSuiteIterator(ParticleSwarmSuite& suite, int updateManager, int topology, int synchronicity);
				bool operator != (ParticleSwarmSuiteIterator const other);
				void operator++(int unused);
				void operator++();
				ParticleSwarm operator * ();
		};

		ParticleSwarmSuiteIterator begin();
		ParticleSwarmSuiteIterator end();
		ParticleSwarm getParticleSwarmByIndex(int const updateManagerIndex, int const topologyIndex, int const synchronicity);

		void setUpdateManagers(std::vector<UpdateManagerType> updateManagers);
		void setTopologyManagers(std::vector<Topology> topologyManagers);
		void setSynchronicities(std::vector<Synchronicity> synchronicities);

		void removeTopology(Topology t);
		void removeUpdateManager(UpdateManagerType u);
		ParticleSwarm getParticleSwarm(int const i);
		std::vector<UpdateManagerType>& getUpdateManagers();
		std::vector<Topology>& getTopologyManagers();
		std::vector<Synchronicity>& getSynchronicities();
		int size() const;
};