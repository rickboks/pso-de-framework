#include "desuite.h"
#include <iostream>
#include <algorithm>

DESuite::DESuite(){
	for (int i = 0; i < MUT_END; i++)
		mutationManagers.push_back((MutationType)i);
	for (int i = 0; i < CROSS_END; i++)
		crossoverManagers.push_back((CrossoverType)i);
	for (int i = 0; i < DEA_END; i++)
		adaptationManagers.push_back((DEAdaptationType)i);
	for (int i = 0; i < INIT_END; i++)
		initializationManagers.push_back((InitializationType)i);

	for (auto mutation : mutationManagers)
		for (auto crossover : crossoverManagers)
				for (auto adaptation : adaptationManagers)
					for (auto initialization : initializationManagers){
						configurations.push_back(std::make_tuple(initialization, mutation, crossover, adaptation, false));
						configurations.push_back(std::make_tuple(initialization, mutation, crossover, adaptation, true));
					}
}


DifferentialEvolution DESuite::getDE(int const i) {
	auto [initialization, mutation, crossover, adapt, jump] = configurations[i];

	return DifferentialEvolution (initialization, mutation, crossover, adapt, jump);
}

void DESuite::setMutationManagers(std::vector<MutationType> mutationManagers){
	this->mutationManagers = mutationManagers;
}

void DESuite::setCrossoverManagers(std::vector<CrossoverType> crossoverManagers){
	this->crossoverManagers = crossoverManagers;
}

void DESuite::setDEAdaptationManagers(std::vector<DEAdaptationType> adaptationManagers){
	this->adaptationManagers = adaptationManagers;
}

void DESuite::setInitializationManagers(std::vector<InitializationType> initialziationManagers){
	this->initializationManagers = initialziationManagers;
}

int DESuite::size() const {
	return configurations.size();
}