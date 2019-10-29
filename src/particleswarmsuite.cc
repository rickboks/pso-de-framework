#include "particleswarmsuite.h"
#include <iostream>
#include <algorithm>

ParticleSwarmSuite::ParticleSwarmSuite(){
	for (int i = 0; i < MAN_END; i++)
		updateManagers.push_back((UpdateManagerType)i);
	for (int i = 0; i < TOP_END; i++)
		topologyManagers.push_back((Topology) i);
	for (int i = 0; i < SYNC_END; i++)
		synchronicities.push_back((Synchronicity)i);

	generateConfigurations();
}

void ParticleSwarmSuite::generateConfigurations(){
	configurations.clear();
	for (auto update : updateManagers)
		for (auto topology : topologyManagers)
			for (auto synchronicity : synchronicities)
					configurations.push_back(
						std::make_tuple(update, topology, synchronicity));

}


ParticleSwarm ParticleSwarmSuite::getParticleSwarm(int const i) {
	auto [update, topology, sync] = configurations[i];

	return ParticleSwarm(update, topology, sync);
}

void ParticleSwarmSuite::setUpdateManagers(std::vector<UpdateManagerType> updateManagers){
	this->updateManagers = updateManagers;
	generateConfigurations();
}

void ParticleSwarmSuite::setTopologyManagers(std::vector<Topology> topologyManagers){
	this->topologyManagers = topologyManagers;
	generateConfigurations();
}

void ParticleSwarmSuite::setSynchronicities(std::vector<Synchronicity> synchronicities){
	this->synchronicities = synchronicities;
	generateConfigurations();
}

int ParticleSwarmSuite::size() const {
	return configurations.size();
}