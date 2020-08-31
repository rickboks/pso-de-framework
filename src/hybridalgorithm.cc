#include "hybridalgorithm.h"

HybridAlgorithm::HybridAlgorithm(UpdateManagerType const updateManagerType, 
			Topology topologyManagerType, Synchronicity const synchronicity, MutationType const mutationType, 
			CrossoverType const crossoverType, SelectionType selection,DEAdaptationType adaptionType):
	updateManagerType(updateManagerType),topologyManagerType(topologyManagerType), 
	synchronicity(synchronicity), mutationType(mutationType), crossoverType(crossoverType), 
	selectionType(selection), adaptationType(adaptionType), topologyManager(NULL), mutationManager(NULL), crossoverManager(NULL),
	adaptationManager(NULL){
}

HybridAlgorithm::~HybridAlgorithm(){

}
