#include "desuite.h"
#include "crossovermanager.h"
#include "deadaptationmanager.h"
#include <iostream>
#include <algorithm>
#include <pthread.h>

DESuite::DESuite(){
	for (auto&i : ::mutations)
		this->mutationManagers.push_back(i.first);
	for (auto&i : ::crossovers)
		this->crossoverManagers.push_back(i.first);
	for (auto&i : ::deAdaptations)
		this->adaptationManagers.push_back(i.first);
	for (auto&i : ::deCHs)
		this->constraintHandlers.push_back(i.first);

	generateConfigurations();
}

void DESuite::generateConfigurations(){
	configurations.clear();
	for (auto mutation : mutationManagers)
		for (auto crossover : crossoverManagers)
				for (auto adaptation : adaptationManagers)
					for (auto ch : constraintHandlers)
							configurations.push_back(DEConfig(mutation, crossover, adaptation, ch));
}

DifferentialEvolution DESuite::getDE(int const i) {
	return DifferentialEvolution (configurations[i]);
}

void DESuite::setMutationManagers(std::vector<std::string> mutationManagers){
	this->mutationManagers = mutationManagers;
	generateConfigurations();
}

void DESuite::setCrossoverManagers(std::vector<std::string> crossoverManagers){
	this->crossoverManagers = crossoverManagers;
	generateConfigurations();
}

void DESuite::setDEAdaptationManagers(std::vector<std::string> adaptationManagers){
	this->adaptationManagers = adaptationManagers;
	generateConfigurations();
}

void DESuite::setConstraintHandlers(std::vector<std::string> chs){
	this->constraintHandlers = chs;
	generateConfigurations();
}

int DESuite::size() const {
	return configurations.size();
}
