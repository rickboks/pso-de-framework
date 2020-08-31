#pragma once
#include <vector>
#include <map>
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "particleupdatemanager.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "hybridalgorithm.h"

template <typename T>
class HybridSuite {
	private:		
		std::vector<UpdateManagerType> updateManagers;
		std::vector<Topology> topologyManagers;
		std::vector<Synchronicity> synchronicities;
		std::vector<MutationType> mutationManagers;
		std::vector<CrossoverType> crossoverManagers;
		std::vector<SelectionType> selectionManagers;
		std::vector<DEAdaptationType> adaptationManagers;
		std::vector<hybrid_config> configurations;

		public:
		HybridSuite(){
			for (int i = 0; i < MUT_END; i++)
				mutationManagers.push_back((MutationType)i);
			for (int i = 0; i < CROSS_END; i++)
				crossoverManagers.push_back((CrossoverType)i);
			for (int i = 0; i < MAN_END; i++)
				updateManagers.push_back((UpdateManagerType)i);
			for (int i = 0; i < TOP_END; i++)
				topologyManagers.push_back((Topology) i);
			for (int i = 0; i < SYNC_END; i++)
				synchronicities.push_back((Synchronicity)i);
			for (int i = 0; i < SEL_END; i++)
				selectionManagers.push_back((SelectionType)i);
			for (int i = 0; i < DEA_END; i++)
				adaptationManagers.push_back((DEAdaptationType)i);

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
										configurations.push_back(
											std::make_tuple(update, topology, synchronicity, mutation, 
												crossover, selection, adaptation));

		}

		T getHybrid(int const i) {
			UpdateManagerType update = std::get<0>(configurations[i]);
			Topology topology = std::get<1>(configurations[i]);
			Synchronicity sync = std::get<2>(configurations[i]);
			MutationType mutation = std::get<3>(configurations[i]);
			CrossoverType crossover = std::get<4>(configurations[i]);
			SelectionType selection = std::get<5> (configurations[i]);
			DEAdaptationType adapt = std::get<6>(configurations[i]);

			return T(update, topology, sync, mutation, crossover, selection, adapt);
		}

		void setMutationManagers(std::vector<MutationType> mutationManagers){
			this->mutationManagers = mutationManagers;
			generateConfigurations();
		}

		void setCrossoverManagers(std::vector<CrossoverType> crossoverManagers){
			this->crossoverManagers = crossoverManagers;
			generateConfigurations();
		}

		void setUpdateManagers(std::vector<UpdateManagerType> updateManagers){
			this->updateManagers = updateManagers;
			generateConfigurations();
		}

		void setTopologyManagers(std::vector<Topology> topologyManagers){
			this->topologyManagers = topologyManagers;
			generateConfigurations();
		}

		void setSynchronicities(std::vector<Synchronicity> synchronicities){
			this->synchronicities = synchronicities;
			generateConfigurations();
		}

		void setSelectionManagers(std::vector<SelectionType> selectionManagers){
			this->selectionManagers = selectionManagers;
			generateConfigurations();
		}

		void setDEAdaptationManagers(std::vector<DEAdaptationType> adaptationManagers){
			this->adaptationManagers = adaptationManagers;
			generateConfigurations();
		}

		int size() const {
			return configurations.size();
		}
};
