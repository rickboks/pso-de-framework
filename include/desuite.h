#pragma once

#include <vector>
#include <map>
#include "differentialevolution.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "deadaptationmanager.h"

typedef std::tuple<DEInitializationType, MutationType, CrossoverType, DEAdaptationType, bool> de_configuration;

class DESuite {
	private:		
		std::vector<MutationType> mutationManagers;
		std::vector<CrossoverType> crossoverManagers;
		std::vector<DEAdaptationType> adaptationManagers;
		std::vector<DEInitializationType> initializationManagers;
		std::vector<de_configuration> configurations;
	public:
		DESuite();
		void generateConfigurations();
		void setMutationManagers(std::vector<MutationType> mutationManagers);
		void setCrossoverManagers(std::vector<CrossoverType> crossoverManagers);
		void setInitializationManagers(std::vector<DEInitializationType> initializationManagers);
		void setDEAdaptationManagers(std::vector<DEAdaptationType> adaptationManagers);
		DifferentialEvolution getDE(int const i);	
		int size() const;
};