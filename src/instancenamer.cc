#include "instancenamer.h"

//For Hybrid
std::string InstanceNamer::getName(UpdateManagerType const update, Topology const topology, 
			Synchronicity const sync, MutationType const mutation, CrossoverType const crossover, 
			SelectionType const selection, DEAdaptationType const adaptation){
	return "H" + updateID(update) + topologyID(topology) + synchronicityID(sync) + mutationID(mutation)
			+ crossoverID(crossover) + selectionID(selection) + adaptationID(adaptation); 
}

//For PSO
std::string InstanceNamer::getName(UpdateManagerType const update, Topology const topology, 
			Synchronicity const sync){
	return "P" + updateID(update) + topologyID(topology) + synchronicityID(sync);
}

//For DE
std::string InstanceNamer::getName(DEInitializationType const init, MutationType const mutation, 
	CrossoverType const crossover, DEAdaptationType const adaptation, bool jump){

	return "D" + deInitializationID(init) + mutationID(mutation) + crossoverID(crossover) + adaptationID(adaptation) + (jump?"1":"0");
}

std::string InstanceNamer::updateID(UpdateManagerType const update){
	switch (update){
		case UpdateManagerType::INERTIA_WEIGHT: return "I";
		case UpdateManagerType::DECR_INERTIA_WEIGHT: return "D";
		case UpdateManagerType::VMAX: return "V";
		case UpdateManagerType::CONSTRICTION_COEFFICIENT: return "C";
		case UpdateManagerType::FIPS: return "F";
		case UpdateManagerType::BARE_BONES: return "B";
		default: return "ERR";
	}
}

std::string InstanceNamer::topologyID(Topology const topology){
	switch (topology){
		case Topology::LBEST: return  "L";
		case Topology::GBEST: return  "G";
		case Topology::RANDOM_GRAPH: return  "R";
		case Topology::VON_NEUMANN:	return  "N";
		case Topology::WHEEL: return  "W";
		case Topology::INCREASING: return  "I";
		case Topology::DECREASING: return  "D";
		case Topology::MULTI_SWARM: return  "M";
		default: return "ERR";
	};	
}

std::string InstanceNamer::synchronicityID(Synchronicity const synchronicity){
	return (synchronicity == SYNCHRONOUS ? "S" : "A");
}

std::string InstanceNamer::mutationID(MutationType const mutation){
	switch (mutation){
		case MutationType::RAND_1: return  "R1";
		case MutationType::BEST_1: return  "B1";
		case MutationType::TTB_1: return  "T1";
		case MutationType::BEST_2: return  "B2";
		case MutationType::RAND_2: return  "R2";
		case MutationType::RAND_2_DIR: return  "RD";
		case MutationType::NSDE: return  "NS";
		case MutationType::TOPOLOGY: return  "TOP";
		case MutationType::TRIGONOMETRIC: return  "TR";
		case MutationType::TTPB_1: return  "PB";
		case MutationType::TO1: return  "O1";
		case MutationType::TO2: return  "O2";
		default: return  "ERR";
	}
	
}

std::string InstanceNamer::crossoverID(CrossoverType const crossover){
	return (crossover == BINOMIAL ? "B" : "E");
}

std::string InstanceNamer::selectionID(SelectionType const selection){
	switch(selection){
		case SelectionType::P2: return  "P2";
		case SelectionType::P3: return  "P3";
		case SelectionType::U2: return  "U2";
		case SelectionType::U3: return  "U3";
		default: return  "ERR";
	}	
}

std::string InstanceNamer::adaptationID(DEAdaptationType const adaptation){
	return (adaptation == JADE ? "J" : "N");
}

std::string InstanceNamer::deInitializationID(DEInitializationType const deinit){
	return (deinit == RANDOM ? "R" : "O");
}