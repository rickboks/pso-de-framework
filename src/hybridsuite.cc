#include "hybridsuite.h"
#include <iostream>
#include <algorithm>

HybridSuite::HybridSuite(){
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

void HybridSuite::generateConfigurations(){
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

HybridAlgorithm HybridSuite::getHybrid(int const i) {
	auto [update, topology, sync, mutation, crossover, selection, adapt] = configurations[i];

	return HybridAlgorithm
	(update, topology, sync, mutation, crossover, selection, adapt);
}

void HybridSuite::setMutationManagers(std::vector<MutationType> mutationManagers){
	this->mutationManagers = mutationManagers;
	generateConfigurations();
}

void HybridSuite::setCrossoverManagers(std::vector<CrossoverType> crossoverManagers){
	this->crossoverManagers = crossoverManagers;
	generateConfigurations();
}

void HybridSuite::setUpdateManagers(std::vector<UpdateManagerType> updateManagers){
	this->updateManagers = updateManagers;
	generateConfigurations();
}

void HybridSuite::setTopologyManagers(std::vector<Topology> topologyManagers){
	this->topologyManagers = topologyManagers;
	generateConfigurations();
}

void HybridSuite::setSynchronicities(std::vector<Synchronicity> synchronicities){
	this->synchronicities = synchronicities;
	generateConfigurations();
}

void HybridSuite::setSelectionManagers(std::vector<SelectionType> selectionManagers){
	this->selectionManagers = selectionManagers;
	generateConfigurations();
}

void HybridSuite::setDEAdaptationManagers(std::vector<DEAdaptationType> adaptationManagers){
	this->adaptationManagers = adaptationManagers;
	generateConfigurations();
}

int HybridSuite::size() const {
	return configurations.size();
}