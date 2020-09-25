#pragma once
#include <vector>
#include <map>
#include "deadaptationmanager.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "particleupdatemanager.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "hybridalgorithm.h"

template <typename T>
class HybridSuite {
	private:		
		std::vector<std::string> updateManagers;
		std::vector<std::string> topologyManagers;
		std::vector<std::string> psoCHs;
		std::vector<std::string> synchronicities;
		std::vector<std::string> mutationManagers;
		std::vector<std::string> crossoverManagers;
		std::vector<std::string> selectionManagers;
		std::vector<std::string> adaptationManagers;
		std::vector<std::string> deCHs;
		std::vector<HybridConfig> configurations;

		public:
		HybridSuite(){
			for (auto&i : ::updateManagers)
				updateManagers.push_back(i.first);

			for (auto&i : ::topologies)
				topologyManagers.push_back(i.first);

			for (auto&i : ::psoCHs)
				psoCHs.push_back(i.first);

			for (auto&i : ::deCHs)
				deCHs.push_back(i.first);

			synchronicities.push_back("A");
			synchronicities.push_back("S");

			for (auto&i : ::mutations)
				mutationManagers.push_back(i.first);
			for (auto&i : ::crossovers)
				crossoverManagers.push_back(i.first);
			//for (auto&i : ::selections)
				//selectionManagers.push_back(i.first);
			for (auto&i : ::deAdaptations)
				adaptationManagers.push_back(i.first);

			generateConfigurations();
		}

		void generateConfigurations(){
			configurations.clear();
			for (auto mutation : mutationManagers)
				for (auto crossover : crossoverManagers)
					for (auto update : updateManagers)
						for (auto topology : topologyManagers)
							for (auto synchronicity : synchronicities)
								for (auto selection : selectionManagers)
									for (auto adaptation : adaptationManagers)
										for (auto deCH : deCHs)
											for (auto psoCH : psoCHs)
												configurations.push_back(HybridConfig(update, topology, psoCH, synchronicity, mutation, 
														crossover, selection, adaptation, deCH));

		}

		T getHybrid(int const i) {
			return T(configurations[i]);
		}

		void setMutationManagers(std::vector<std::string> mutationManagers){
			this->mutationManagers = mutationManagers;
			generateConfigurations();
		}

		void setCrossoverManagers(std::vector<std::string> crossoverManagers){
			this->crossoverManagers = crossoverManagers;
			generateConfigurations();
		}

		void setUpdateManagers(std::vector<std::string> updateManagers){
			this->updateManagers = updateManagers;
			generateConfigurations();
		}

		void setTopologyManagers(std::vector<std::string> topologyManagers){
			this->topologyManagers = topologyManagers;
			generateConfigurations();
		}

		void setSynchronicities(std::vector<std::string> synchronicities){
			this->synchronicities = synchronicities;
			generateConfigurations();
		}

		void setSelectionManagers(std::vector<std::string> selectionManagers){
			this->selectionManagers = selectionManagers;
			generateConfigurations();
		}

		void setDEAdaptationManagers(std::vector<std::string> adaptationManagers){
			this->adaptationManagers = adaptationManagers;
			generateConfigurations();
		}
		void setDEContraintHandlers(std::vector<std::string> deCHs){
			this->deCHs = deCHs;
			generateConfigurations();
		}

		void setPSOConstraintHandlers(std::vector<std::string> psoCHs){
			this->psoCHs = psoCHs;
			generateConfigurations();
		}

		int size() const {
			return configurations.size();
		}
};
