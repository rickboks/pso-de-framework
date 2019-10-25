#include "particleupdatemanager.h"
#include "topologymanager.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "deadaptationmanager.h"
#include "selectionmanager.h"
#include "particleswarm.h"
#include "deinitializer.h"

class InstanceNamer {
	private:
		static std::string updateID(UpdateManagerType const update);
		static std::string topologyID(Topology const topology);
		static std::string synchronicityID(Synchronicity const synchronicity);
		static std::string mutationID(MutationType const mutation);
		static std::string crossoverID(CrossoverType const crossover);
		static std::string selectionID(SelectionType const selection);
		static std::string adaptationID(DEAdaptationType const adaptation);
		static std::string deInitializationID(DEInitializationType deinit);
	public:
		//For hyrbid
		static std::string getName(UpdateManagerType const update, Topology const topology, 
			Synchronicity const sync, MutationType const mutation, CrossoverType const crossover, 
			SelectionType const selection, DEAdaptationType const adaptation);

		//For PSO
		static std::string getName(UpdateManagerType const update, Topology const topology, 
			Synchronicity const sync);

		//For DE
		static std::string getName(DEInitializationType const init, MutationType const mutation, CrossoverType const crossover, 
			DEAdaptationType const adaptation, bool jump);

};