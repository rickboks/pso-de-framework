#include "particleswarmsuite.h"
#include <iostream>
#include <algorithm>
ParticleSwarmSuite::ParticleSwarmSuite(std::map<int,double>& updateParams)
	:updateParams(updateParams){

	for (int i = 0; i < TOP_END; i++){
		topologyManagers.push_back((Topology)i);
	}

	for (int i = 0; i < MAN_END; i++){
		updateManagers.push_back((UpdateManagerType)i);
	}

	for (int i = 0; i < SYNC_END; i++){
		synchronicities.push_back((Synchronicity)i);
	}
}

ParticleSwarmSuite::ParticleSwarmSuiteIterator ParticleSwarmSuite::begin(){
	return ParticleSwarmSuite::ParticleSwarmSuiteIterator(*this, 0, 0, 0);
}

ParticleSwarmSuite::ParticleSwarmSuiteIterator ParticleSwarmSuite::end(){
	return ParticleSwarmSuite::ParticleSwarmSuiteIterator(*this, updateManagers.size(), topologyManagers.size(), synchronicities.size());
}

ParticleSwarm ParticleSwarmSuite::getParticleSwarmByIndex(int const updateManagerIndex, 
	int const topologyIndex, int const synchronicity){

	return ParticleSwarm(
		updateManagers[updateManagerIndex], 
		topologyManagers[topologyIndex], synchronicities[synchronicity]);
}

ParticleSwarm ParticleSwarmSuite::getParticleSwarm(int const i){
	ParticleSwarmSuite::ParticleSwarmSuiteIterator it = begin();

	for (int n = 0; n < i; n++)
		it++;

	return *it;
}

void ParticleSwarmSuite::setUpdateManagers(std::vector<UpdateManagerType> updateManagers){
	this->updateManagers = updateManagers;
}

void ParticleSwarmSuite::setTopologyManagers(std::vector<Topology> topologyManagers){
	this->topologyManagers = topologyManagers;
}

void ParticleSwarmSuite::setSynchronicities(std::vector<Synchronicity> synchronicities){
	this->synchronicities = synchronicities;
}

std::vector<UpdateManagerType>& ParticleSwarmSuite::getUpdateManagers(){
	return updateManagers;
}

std::vector<Topology>& ParticleSwarmSuite::getTopologyManagers(){
	return topologyManagers;
}

std::vector<Synchronicity>& ParticleSwarmSuite::getSynchronicities(){
	return synchronicities;
}

/// ITERATOR
ParticleSwarmSuite::ParticleSwarmSuiteIterator::ParticleSwarmSuiteIterator( 
	ParticleSwarmSuite& suite, int updateManager, int topology, int synchronicity)
	:updateManagers(suite.getUpdateManagers()), topologyManagers(suite.getTopologyManagers()),
	synchronicities(suite.getSynchronicities()), updateManager(updateManager), topology(topology), synchronicity(synchronicity){
}

bool ParticleSwarmSuite::ParticleSwarmSuiteIterator::operator != (ParticleSwarmSuiteIterator const other){
	return (other.updateManager != this->updateManager ||
			other.topology != this->topology ||
			other.synchronicity != this->synchronicity);
} 

void ParticleSwarmSuite::ParticleSwarmSuiteIterator::operator++(int unused){
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
				synchronicity = synchronicities.size();
				topology = topologyManagers.size();
				updateManager = updateManagers.size();
			}
		}
	}
}

void ParticleSwarmSuite::ParticleSwarmSuiteIterator::operator++(){
	++(*this);
}

ParticleSwarm ParticleSwarmSuite::ParticleSwarmSuiteIterator::operator * (){
		return ParticleSwarm(
		updateManagers[updateManager], 
		topologyManagers[topology], synchronicities[synchronicity]);
}

int ParticleSwarmSuite::size() const {
	return updateManagers.size() * topologyManagers.size() * synchronicities.size();
}


void ParticleSwarmSuite::removeTopology(Topology t){
	topologyManagers.erase(std::remove(topologyManagers.begin(), topologyManagers.end(), t), topologyManagers.end());
}
void ParticleSwarmSuite::removeUpdateManager(UpdateManagerType u){
	updateManagers.erase(std::remove(updateManagers.begin(), updateManagers.end(), u), updateManagers.end());
}