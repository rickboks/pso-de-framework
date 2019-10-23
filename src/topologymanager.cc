#include "topologymanager.h"
#include "particle.h"
#include "rng.h"
#include <algorithm>
#include <iostream>

/*		Base 		*/
TopologyManager::TopologyManager(std::vector<Particle*> const & particles)
	:particles(particles){	
}

TopologyManager::~TopologyManager() = default;

int TopologyManager::getClosestValidPopulationSize(int const popSize){
	if (popSize > 0)
		return popSize;
	else {
		std::cerr << "Warning: population size needs to be greater or equal to 1. Setting population size to 1." << std::endl;
		return 1;
	}
}

void TopologyManager::update(double progress){
	// default: do nothing. (static topologies)
}


/*		Factory 	*/
TopologyManager* TopologyManagerFactory::createTopologyManager(int const type, std::vector<Particle*> const & particles){
	switch(type){
		case Topology::LBEST:
			return new LbestTopologyManager(particles);
		case Topology::GBEST:
			return new GbestTopologyManager(particles);
		case Topology::RANDOM_GRAPH:
			return new RandomTopologyManager(particles);
		case Topology::VON_NEUMANN:
			return new VonNeumannTopologyManager(particles);
		case Topology::WHEEL:
			return new WheelTopologyManager(particles);
		case Topology::INCREASING:
			return new IncreasingTopologyManager(particles);
		case Topology::DECREASING:
			return new DecreasingTopologyManager(particles);
		case Topology::MULTI_SWARM:
			return new MultiSwarmTopologyManager(particles);
		default:
			throw std::invalid_argument("Error: Invalid topology manager type");
	}
}

/*		Lbest 		*/
LbestTopologyManager::LbestTopologyManager(std::vector<Particle*> const & particles)
	:TopologyManager(particles){

}

void LbestTopologyManager::initialize(){
	int const popSize = particles.size();
	for (int i = 0; i < popSize; i++){
		if (i ==0)
			particles[0]->addNeighbor(particles[popSize-1]);
		else
			particles[i]->addNeighbor(particles[i-1]);

		if (i == popSize -1)
			particles[i]->addNeighbor(particles[0]);
		else
			particles[i]->addNeighbor(particles[i+1]);
	}
}

int LbestTopologyManager::getClosestValidPopulationSize(int const popSize){
	if (popSize >= 3)
		return popSize;
	else {
		std::cerr << "Warning: Lbest topology needs a population size greater or equal to " << 3 <<
		". Setting population size to " << 3 << "." << std::endl;
		return 3;
	}
}

/*		Gbest 		*/
GbestTopologyManager::GbestTopologyManager(std::vector<Particle*> const & particles)
	:TopologyManager(particles){

}

void GbestTopologyManager::initialize(){
	int const popSize = particles.size();
	for (int i = 0; i < popSize; i++){
		for (int j = 0; j < popSize; j++){
			if (i != j){
				particles[i]->addNeighbor(particles[j]);
			}
		}
	}
}

/*		Random 		*/
RandomTopologyManager::RandomTopologyManager(std::vector<Particle*> const & particles)
	:TopologyManager(particles), connections(3){

}

void RandomTopologyManager::initialize(){
	int const popSize = particles.size();

	for (int i = 0; i < popSize; i++){
		std::vector<int> possibilities;
		possibilities.reserve(popSize-1);
		for (int j = 0; j < popSize; j++){
			if (i != j){
				possibilities.push_back(j);
			}
		}

		for (int j = 0; j < connections; j++){
			int index = rng.randInt(0,possibilities.size()-1);
			particles[i]->addNeighbor(particles[possibilities[index]]);
			possibilities.erase(possibilities.begin() + index);
		}
	}
}

int RandomTopologyManager::getClosestValidPopulationSize(int const popSize){
	if (popSize > connections){
		return popSize;
	} else {
		std::cerr << "Warning: Random graph topology needs a population size of at least " << connections+1 << ". "
		<< "Setting population size to " << connections +1 << "." << std::endl;
		return connections+1;
	}
}

/*		Von Neumann		*/
VonNeumannTopologyManager::VonNeumannTopologyManager(std::vector<Particle*> const & particles)
	:TopologyManager(particles){
}

void VonNeumannTopologyManager::initialize(){
	int rows = (int)sqrt(particles.size());
	  while (particles.size() % rows != 0) {
	    rows--;
  	}
  	int const columns = particles.size() / rows;

	for (int i = 0; i < rows; i++){
		for (int j = 0; j < columns; j++){
			if (i != rows-1)
				particles[i * columns + j]->addNeighbor(particles[(i+1) * columns + j]);
			else 
				particles[i * columns + j]->addNeighbor(particles[j]);

			if (i != 0)
				particles[i * columns + j]->addNeighbor(particles[(i-1) * columns + j]);
			else 
				particles[i * columns + j]->addNeighbor(particles[(rows-1)* columns + j]);

			if (j != columns-1)
				particles[i * columns + j]->addNeighbor(particles[i * columns + (j+1)]);
			else 
				particles[i * columns + j]->addNeighbor(particles[i * columns + 0]);

			if (j != 0)
				particles[i * columns + j]->addNeighbor(particles[i * columns + (j-1)]);
			else 
				particles[i * columns + j]->addNeighbor(particles[i * columns + (columns -1)]);
		}
	}
}

/*		Wheel 		*/
WheelTopologyManager::WheelTopologyManager(std::vector<Particle*> const & particles)
	:TopologyManager(particles){
}

void WheelTopologyManager::initialize(){
	int const popSize = particles.size();
	for (int i = 1; i < popSize; i++){
		particles[i]->addNeighbor(particles[0]);
		particles[0]->addNeighbor(particles[i]);
	}
}

/*		Increasing connectivity 	*/
IncreasingTopologyManager::IncreasingTopologyManager(std::vector<Particle*> const & particles)
	:TopologyManager(particles), currentConnectivity(2), minConnectivity(2){
}

void IncreasingTopologyManager::initialize(){
	maxConnectivity = particles.size() -1;

	int const popSize = particles.size();
	for (int i = 0; i < popSize; i++){
		if (i ==0)
			particles[0]->addNeighbor(particles[popSize-1]);
		else
			particles[i]->addNeighbor(particles[i-1]);

		if (i == popSize -1)
			particles[i]->addNeighbor(particles[0]);
		else
			particles[i]->addNeighbor(particles[i+1]);
	}


}

void IncreasingTopologyManager::update(double progress){
	progress = std::min(1.0, progress);
	int const popSize = particles.size();
	int const newConnectivity = minConnectivity + progress * (maxConnectivity - minConnectivity);

	if (newConnectivity > currentConnectivity){
		int const newNeighbors = newConnectivity - currentConnectivity;
		std::vector<Particle*> possibilities;

		for (int i = 0; i < popSize; i++){
			for (int k = 0; k < popSize; k++){
				if (k != i && !particles[i]->isNeighbor(particles[k])){
					possibilities.push_back(particles[k]);
				}
			}

			for (int j = 0; j < newNeighbors; j++){
				int randIndex = rng.randInt(0,possibilities.size()-1);
				particles[i]->addNeighbor(possibilities[randIndex]);
				possibilities.erase(possibilities.begin() + randIndex);			
			}			
			possibilities.clear();
			
		}
		currentConnectivity = newConnectivity;
	}
}

/* Decreasing connectivity */

DecreasingTopologyManager::DecreasingTopologyManager(std::vector<Particle*> const & particles)
	:TopologyManager(particles){
}

void DecreasingTopologyManager::initialize(){
	int const popSize = particles.size();
	maxConnectivity = popSize;
	currentConnectivity = popSize-1;
	for (int i = 0; i < popSize; i++){
		for (int j = 0; j < popSize; j++){
			if (i != j)
				particles[i]->addNeighbor(particles[j]);
		}
	}
}

void DecreasingTopologyManager::update(double progress){
	progress = std::min(1.0, progress);
	int const popSize = particles.size();
	int const newConnectivity = maxConnectivity - progress * (maxConnectivity - 2);

	if (newConnectivity < currentConnectivity){
		int const removeNeighbors = currentConnectivity - newConnectivity;
		std::vector<Particle*> possibilities;

		for (int i = 0; i < popSize; i++){
			for (int k = 0; k < popSize; k++){
				if (i != k && !(i == 0 && k == popSize-1) && k != i-1  && k != ((i+1)%popSize)
					&& particles[i]->isNeighbor(particles[k])){
					possibilities.push_back(particles[k]);
				}
			}

			for (int k = 0; k < removeNeighbors; k++){
				int randIndex = rng.randInt(0,possibilities.size()-1);
				particles[i]->removeNeighbor(possibilities[randIndex]);
				possibilities.erase(possibilities.begin() + randIndex);

			}
			possibilities.clear();
			
		}
		currentConnectivity = newConnectivity;
	}
}


/* Dynamic multi-swarm */

MultiSwarmTopologyManager::MultiSwarmTopologyManager(std::vector<Particle*> const & particles)
	:TopologyManager(particles), clusterSize(3), count(0){

}

void MultiSwarmTopologyManager::createClusters(){
	std::vector<Particle*> toInitialize = particles;
	
	while (!toInitialize.empty()){
		std::vector<Particle*> cluster;
		int newClusterSize;

		if ((int)toInitialize.size() >= 2 * clusterSize)
			newClusterSize = clusterSize;
		else 
			newClusterSize = toInitialize.size();

		cluster.reserve(newClusterSize);
		for (int i = 0; i < newClusterSize; i++){
			int randIndex = rng.randInt(0,toInitialize.size()-1);
			cluster.push_back(toInitialize[randIndex]);
			toInitialize.erase(toInitialize.begin() + randIndex);
		}

		for (int i = 0; i < (int)cluster.size(); i++){
			for (int j = 0; j < (int)cluster.size(); j++){
				if (i != j)
					cluster[i]->addNeighbor(cluster[j]);
			}
		}

		cluster.clear();
	}
}

void MultiSwarmTopologyManager::initialize(){
	createClusters();
}

void MultiSwarmTopologyManager::update(double progress){
	count++;
	if (count >= 5){
		for (Particle* p: particles)
			p->removeAllNeighbors();
		createClusters();
		count = 0;
	}
}

int MultiSwarmTopologyManager::getClosestValidPopulationSize(int const popSize){
	if (popSize >= clusterSize){
		return popSize;
	} else {
		std::cerr << "Warning: multiswarm topology needs a population size of at least " << clusterSize << ". "
		<< "Setting population size to " << clusterSize << "." << std::endl;
		return clusterSize;
	}
}
