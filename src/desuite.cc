#include "desuite.h"
#include <iostream>
#include <algorithm>

DESuite::DESuite(){
	for (int i = 0; i < INIT_END; i++){
		initializationManagers.push_back((InitializationType)i);
	}

	for (int i = 0; i < MUT_END; i++){
		mutationManagers.push_back((MutationType)i);
	}

	for (int i = 0; i < CROSS_END; i++){
		crossoverManagers.push_back((CrossoverType)i);
	}
}

DESuite::DESuiteIterator DESuite::begin(){
	return DESuite::DESuiteIterator(*this, 0, 0, 0, false);
}

DESuite::DESuiteIterator DESuite::end(){
	return DESuite::DESuiteIterator(*this, initializationManagers.size(), mutationManagers.size(), crossoverManagers.size(), true);
}

DifferentialEvolution DESuite::getDEByIndex(int const intializationManagerIndex, 
	int const mutationManagerIndex, int const crossoverManagerIndex, bool const jumpOpposition) const{

	return DifferentialEvolution(
		initializationManagers[intializationManagerIndex], 
		mutationManagers[mutationManagerIndex], crossoverManagers[crossoverManagerIndex], jumpOpposition);
}

DifferentialEvolution DESuite::getDE(int const i) {
	DESuite::DESuiteIterator it = begin();


	for (int n = 0; n < i; n++){
		it++;
	}
	return *it;
}

void DESuite::setInitializationManagers(std::vector<InitializationType> initializationManagers){
	this->initializationManagers = initializationManagers;
}

void DESuite::setMutationManagers(std::vector<MutationType> mutationManagers){
	this->mutationManagers = mutationManagers;
}

void DESuite::setCrossoverManagers(std::vector<CrossoverType> crossoverManagers){
	this->crossoverManagers = crossoverManagers;
}


std::vector<InitializationType> DESuite::getInitializationManagers() const{
	return initializationManagers;
}

std::vector<MutationType> DESuite::getMutationManagers() const{
	return mutationManagers;
}

std::vector<CrossoverType> DESuite::getCrossoverManagers() const{
	return crossoverManagers;
}

/// ITERATOR
DESuite::DESuiteIterator::DESuiteIterator( 
	DESuite& suite, int initialization, int mutation, int crossover, bool jumpOpposition)
	:initializationManagers(suite.getInitializationManagers()), mutationManagers(suite.getMutationManagers()),
	crossoverManagers(suite.getCrossoverManagers()), initialization(initialization), 
	mutation(mutation), crossover(crossover), jumpOpposition(jumpOpposition){
}

bool DESuite::DESuiteIterator::operator != (DESuiteIterator const other){
	return (other.initialization != this->initialization ||
			other.mutation != this->mutation ||
			other.crossover != this->crossover ||
			other.jumpOpposition != this->jumpOpposition);
} 

void DESuite::DESuiteIterator::operator++(int unused){
	if (initialization != (int)initializationManagers.size() -1){
		initialization++;
	} else {
		initialization = 0;
		if (mutation != (int)mutationManagers.size() -1) {
			mutation++;
		} else {
			mutation = 0;
			if (crossover != (int)crossoverManagers.size()-1) {
				crossover++;
			} else {
				crossover = 0;
				if (jumpOpposition != true){
					jumpOpposition = true;
				} else {
					initialization = initializationManagers.size();
					mutation = mutationManagers.size();
					crossover = crossoverManagers.size();
				}
			}
		}
	}
}

void DESuite::DESuiteIterator::operator++(){
	++(*this);
}

DifferentialEvolution DESuite::DESuiteIterator::operator * () {
	return DifferentialEvolution(
		initializationManagers[initialization], 
		mutationManagers[mutation], crossoverManagers[crossover], jumpOpposition);
}

int DESuite::size() const {
	return initializationManagers.size() * mutationManagers.size() * crossoverManagers.size() * 2;
}