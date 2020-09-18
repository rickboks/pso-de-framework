#include "particleswarmsuite.h"
#include <iostream>
#include <algorithm>

ParticleSwarmSuite::ParticleSwarmSuite(){
	for (auto& i : ::updateManagers)
		this->updateManagers.push_back(i.first);
	for (auto& i : ::topologies)
		this->topologyManagers.push_back(i.first);
	for (auto& i : ::psoCHs)
		this->constraintHandlers.push_back(i.first);

	this->synchronicities.push_back("A");
	this->synchronicities.push_back("S");

	generateConfigurations();
}

void ParticleSwarmSuite::generateConfigurations(){
	configurations.clear();
	for (auto update : updateManagers)
		for (auto topology : topologyManagers)
			for (auto ch : constraintHandlers)
				for (auto synchronicity : synchronicities)
					configurations.push_back(PSOConfig(update, topology, ch, synchronicity));
}

ParticleSwarm ParticleSwarmSuite::getParticleSwarm(int const i) {
	return ParticleSwarm(configurations[i]);
}

void ParticleSwarmSuite::setUpdateManagers(std::vector<std::string> updateManagers){
	this->updateManagers = updateManagers;
	generateConfigurations();
}

void ParticleSwarmSuite::setTopologyManagers(std::vector<std::string> topologyManagers){
	this->topologyManagers = topologyManagers;
	generateConfigurations();
}

void ParticleSwarmSuite::setConstraintHandlers(std::vector<std::string> chs){
	this->constraintHandlers = chs;
	generateConfigurations();
}

void ParticleSwarmSuite::setSynchronicities(std::vector<std::string> synchronicities){
	this->synchronicities = synchronicities;
	generateConfigurations();
}

int ParticleSwarmSuite::size() const {
	return configurations.size();
}
