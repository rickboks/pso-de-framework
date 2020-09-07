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

	generateConfigurations();
}

void DESuite::generateConfigurations(){
	configurations.clear();
	for (auto mutation : mutationManagers)
		for (auto crossover : crossoverManagers)
				for (auto adaptation : adaptationManagers)
						configurations.push_back(std::make_tuple(mutation, crossover, adaptation));
						//configurations.push_back(std::make_tuple(initialization, mutation, crossover, adaptation, true));

}

DifferentialEvolution DESuite::getDE(int const i) {
	//auto [initialization, mutation, crossover, adapt, jump] = configurations[i];
	MutationType mutation = std::get<0>(configurations[i]);
	CrossoverType crossover = std::get<1>(configurations[i]);
	DEAdaptationType adapt = std::get<2>(configurations[i]);
	return DifferentialEvolution (mutation, crossover, adapt);
}

void DESuite::setMutationManagers(std::vector<MutationType> mutationManagers){
	this->mutationManagers = mutationManagers;
	generateConfigurations();
}

void DESuite::setCrossoverManagers(std::vector<CrossoverType> crossoverManagers){
	this->crossoverManagers = crossoverManagers;
	generateConfigurations();
}

void DESuite::setDEAdaptationManagers(std::vector<DEAdaptationType> adaptationManagers){
	this->adaptationManagers = adaptationManagers;
	generateConfigurations();
}

int DESuite::size() const {
	return configurations.size();
}
