#include "hybridsuite.h"
#include <iostream>
#include <algorithm>

HybridSuite::HybridSuite(std::map<int,double> & updateParams)
	:updateParams(updateParams){

	for (int i = 0; i < MUT_END; i++){
		mutationManagers.push_back((MutationType)i);
	}

	for (int i = 0; i < CROSS_END; i++){
		crossoverManagers.push_back((CrossoverType)i);
	}

	for (int i = 0; i < MAN_END; i++){
		updateManagers.push_back((UpdateManagerType)i);
	}

	for (int i = 0; i < TOP_END; i++){
		topologyManagers.push_back((Topology) i);
	}

	for (int i = 0; i < SYNC_END; i++){
		synchronicities.push_back((Synchronicity)i);
	}
}

HybridSuite::HybridSuiteIterator HybridSuite::begin(){
	return HybridSuite::HybridSuiteIterator(*this, 0, 0, 0, 0, 0);
}

HybridSuite::HybridSuiteIterator HybridSuite::end(){
	return HybridSuite::HybridSuiteIterator(*this, updateManagers.size(), topologyManagers.size(), synchronicities.size(), mutationManagers.size(), crossoverManagers.size());
}

HybridAlgorithm HybridSuite::getHybridByIndex(int const update, int const topology, int const synchronicity, 
	int const mutation, int const crossover) const{

	return HybridAlgorithm(updateManagers[update], topologyManagers[topology], synchronicities[synchronicity], 
		mutationManagers[mutation], crossoverManagers[crossover]);
}

HybridAlgorithm HybridSuite::getHybrid(int const i) {
	HybridSuite::HybridSuiteIterator it = begin();


	for (int n = 0; n < i; n++){
		it++;
	}
	return *it;
}

void HybridSuite::setMutationManagers(std::vector<MutationType> mutationManagers){
	this->mutationManagers = mutationManagers;
}

void HybridSuite::setCrossoverManagers(std::vector<CrossoverType> crossoverManagers){
	this->crossoverManagers = crossoverManagers;
}

void HybridSuite::setUpdateManagers(std::vector<UpdateManagerType> updateManagers){
	this->updateManagers = updateManagers;
}

void HybridSuite::setTopologyManagers(std::vector<Topology> topologyManagers){
	this->topologyManagers = topologyManagers;
}

void HybridSuite::setSynchronicities(std::vector<Synchronicity> synchronicities){
	this->synchronicities = synchronicities;
}

std::vector<UpdateManagerType> HybridSuite::getUpdateManagers() const{
	return updateManagers;
}

std::vector<Topology> HybridSuite::getTopologyManagers() const{
	return topologyManagers;
}

std::vector<Synchronicity>  HybridSuite::getSynchronicities() const{
	return synchronicities;
}

std::vector<MutationType> HybridSuite::getMutationManagers() const{
	return mutationManagers;
}

std::vector<CrossoverType> HybridSuite::getCrossoverManagers() const{
	return crossoverManagers;
}


/// ITERATOR
HybridSuite::HybridSuiteIterator::HybridSuiteIterator( 
	HybridSuite& suite, int update, int topology, int synchronicity,  int mutation, int crossover)
	:updateManagers(suite.getUpdateManagers()), topologyManagers(suite.getTopologyManagers()), 
	synchronicities(suite.getSynchronicities()), mutationManagers(suite.getMutationManagers()),
	crossoverManagers(suite.getCrossoverManagers()), updateManager(update), topology(topology),
	synchronicity(synchronicity), mutation(mutation), crossover(crossover) {
}

bool HybridSuite::HybridSuiteIterator::operator != (HybridSuiteIterator const other){
	return (other.updateManager != this->updateManager ||
			other.topology != this->topology ||
			other.synchronicity != this->synchronicity ||
			other.mutation != this->mutation ||
			other.crossover != this->crossover);
} 

void HybridSuite::HybridSuiteIterator::operator++(int unused){

	if (updateManager != (int)updateManagers.size() -1){
		updateManager++;
	} else {
		updateManager = 0;
		if (topology != (int)topologyManagers.size() -1) {
			topology++;
		} else {
			topology = 0;
			if (synchronicity != (int)synchronicities.size()-1) {
				synchronicity++;
			} else {
				synchronicity = 0;
				if (mutation != (int)mutationManagers.size() -1){
					mutation++;
				} else {
					mutation = 0;
					if (crossover != (int) crossoverManagers.size() -1){
						crossover++;
					} else {
						updateManager = updateManagers.size();
						topology = topologyManagers.size();
						synchronicity = synchronicities.size();
						mutation = mutationManagers.size();
						crossover = crossoverManagers.size();
					}					
				}
			}
		}
	}
}

void HybridSuite::HybridSuiteIterator::operator++(){
	if (updateManager != (int)updateManagers.size() -1){
		updateManager++;
	} else {
		updateManager = 0;
		if (topology != (int)topologyManagers.size() -1) {
			topology++;
		} else {
			topology = 0;
			if (synchronicity != (int)synchronicities.size()-1) {
				synchronicity++;
			} else {
				synchronicity = 0;
				if (mutation != (int)mutationManagers.size() -1){
					mutation++;
				} else {
					mutation = 0;
					if (crossover != (int)crossoverManagers.size() -1){
						crossover++;
					} else {
						updateManager = updateManagers.size();
						topology = topologyManagers.size();
						synchronicity = synchronicities.size();
						mutation = mutationManagers.size();
						crossover = crossoverManagers.size();
					}					
				}
			}
		}
	}
}

HybridAlgorithm HybridSuite::HybridSuiteIterator::operator * () {
	return HybridAlgorithm(
		updateManagers[updateManager], topologyManagers[topology], synchronicities[synchronicity],
		mutationManagers[mutation], crossoverManagers[crossover]);
}

int HybridSuite::size() const {
	return updateManagers.size() * topologyManagers.size() * synchronicities.size() * 
		mutationManagers.size() * crossoverManagers.size();
}