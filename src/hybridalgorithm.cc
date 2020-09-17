#include "hybridalgorithm.h"

HybridAlgorithm::HybridAlgorithm(UpdateManagerType const updateManagerType, 
	Topology topologyManagerType, Synchronicity const synchronicity, MutationType const mutationType, 
	CrossoverType const crossoverType, SelectionType const selection,DEAdaptationType const adaptionType, 
	std::string const psoCH, std::string const deCH)
		: config(updateManagerType, topologyManagerType, synchronicity, 
		mutationType, crossoverType, selection, adaptionType, psoCH, deCH){}

HybridAlgorithm::HybridAlgorithm(hybrid_config config)
	: config(config){}

HybridAlgorithm::~HybridAlgorithm(){}
