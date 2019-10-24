#pragma once
#include <vector>
#include <map>
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "particleupdatemanager.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "hybridalgorithm.h"
typedef std::tuple<UpdateManagerType, Topology, Synchronicity, MutationType, 
	CrossoverType, SelectionType, DEAdaptationType> hybrid_configuration;

class HybridSuite {
	private:		
		std::vector<UpdateManagerType> updateManagers;
		std::vector<Topology> topologyManagers;
		std::vector<Synchronicity> synchronicities;
		std::vector<MutationType> mutationManagers;
		std::vector<CrossoverType> crossoverManagers;
		std::vector<SelectionType> selectionManagers;
		std::vector<DEAdaptationType> adaptationManagers;
		std::vector<hybrid_configuration> configurations;
		public:
		HybridSuite();
		void setMutationManagers(std::vector<MutationType> mutationManagers);
		void setCrossoverManagers(std::vector<CrossoverType> crossoverManagers);
		void setUpdateManagers(std::vector<UpdateManagerType> updateManagers);
		void setTopologyManagers(std::vector<Topology> topologyManagers);
		void setSynchronicities(std::vector<Synchronicity> synchronicities);
		void setSelectionManagers(std::vector<SelectionType> selectionManagers);
		void setDEAdaptationManagers(std::vector<DEAdaptationType> adaptationManagers);
		void generateConfigurations();
		HybridAlgorithm getHybrid(int const i);	
		int size() const;
};