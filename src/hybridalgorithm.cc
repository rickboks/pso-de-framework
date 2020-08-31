#include "hybridalgorithm.h"

HybridAlgorithm::HybridAlgorithm(UpdateManagerType const updateManagerType, 
			Topology topologyManagerType, Synchronicity const synchronicity, MutationType const mutationType, 
			CrossoverType const crossoverType, SelectionType selection,DEAdaptationType adaptionType): config(updateManagerType, topologyManagerType, synchronicity, mutationType, crossoverType, selection, adaptionType)
	{
}

HybridAlgorithm::~HybridAlgorithm(){

}
