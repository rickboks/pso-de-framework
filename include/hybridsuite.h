#pragma once

#include <vector>
#include <map>
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "particleupdatemanager.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "hybridalgorithm.h"
class HybridSuite {
	private:		
		std::map<int,double>& updateParams;
		std::vector<UpdateManagerType> updateManagers;
		std::vector<Topology> topologyManagers;
		std::vector<Synchronicity> synchronicities;
		std::vector<MutationType> mutationManagers;
		std::vector<CrossoverType> crossoverManagers;		
	public:
		HybridSuite(std::map<int,double> & updateParams);

		class HybridSuiteIterator {
			private:
				std::vector<UpdateManagerType> updateManagers;
				std::vector<Topology> topologyManagers;
				std::vector<Synchronicity> synchronicities;
				std::vector<MutationType> mutationManagers;
				std::vector<CrossoverType> crossoverManagers;				
				int updateManager;
				int topology;
				int synchronicity;
				int mutation;
				int crossover;
			public:
				HybridSuiteIterator (HybridSuite& suite, int update, int topology, int synchronicity, int mutation, int crossover);
				bool operator != (HybridSuiteIterator const other);
				void operator++(int unused);
				void operator++();
				HybridAlgorithm operator *();
		};

		HybridSuiteIterator begin();
	    HybridSuiteIterator end();
		HybridAlgorithm getHybridByIndex(int const update, int const topology, int const synchronicity, 
			int const mutation, int const crossover) const;

		void setMutationManagers(std::vector<MutationType> mutationManagers);
		void setCrossoverManagers(std::vector<CrossoverType> crossoverManagers);
		void setUpdateManagers(std::vector<UpdateManagerType> updateManagers);
		void setTopologyManagers(std::vector<Topology> topologyManagers);
		void setSynchronicities(std::vector<Synchronicity> synchronicities);
		HybridAlgorithm getHybrid(int const i);
		std::vector<MutationType> getMutationManagers() const;
		std::vector<CrossoverType> getCrossoverManagers() const;
		std::vector<UpdateManagerType> getUpdateManagers() const;
		std::vector<Topology> getTopologyManagers() const;
		std::vector<Synchronicity> getSynchronicities() const;
	
		int size() const;
};