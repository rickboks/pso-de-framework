#pragma once

#include <vector>
#include <map>
#include "differentialevolution.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "deadaptationmanager.h"

class DESuite {
	private:		
		std::vector<std::string> mutationManagers;
		std::vector<std::string> crossoverManagers;
		std::vector<std::string> adaptationManagers;
		std::vector<std::string> constraintHandlers;
		std::vector<DEConfig> configurations;
	public:
		DESuite();
		void generateConfigurations();
		void setMutationManagers(std::vector<std::string> mutationManagers);
		void setCrossoverManagers(std::vector<std::string> crossoverManagers);
		void setDEAdaptationManagers(std::vector<std::string> adaptationManagers);
		void setConstraintHandlers(std::vector<std::string> adaptationManagers);
		DifferentialEvolution getDE(int const i);	
		int size() const;
};
