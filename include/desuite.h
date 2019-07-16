#pragma once

#include <vector>
#include "deinitializer.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "differentialevolution.h"
class DESuite {
	private:		
		std::vector<InitializationType> initializationManagers;
		std::vector<MutationType> mutationManagers;
		std::vector<CrossoverType> crossoverManagers;
		bool jumpOpposition;
	public:
		DESuite();

		class DESuiteIterator {
			private:
				std::vector<InitializationType> const initializationManagers;
				std::vector<MutationType> const mutationManagers;
				std::vector<CrossoverType> const crossoverManagers;
				int initialization;
				int mutation;
				int crossover;
				bool jumpOpposition;
			public:
				DESuiteIterator (DESuite& suite, int initType, int mutation, int crossover, bool const jumpOpposition);
				bool operator != (DESuiteIterator const other);
				void operator++(int unused);
				void operator++();
				DifferentialEvolution operator *();
		};

		DESuiteIterator begin();
	    DESuiteIterator end();
		DifferentialEvolution getDEByIndex(int const initIndex, int const mutationIndexIndex, int const crossoverIndex, bool const jumpOpposition) const;

		void setInitializationManagers(std::vector<InitializationType> initializationManagers);
		void setMutationManagers(std::vector<MutationType> mutationManagers);
		void setCrossoverManagers(std::vector<CrossoverType> crossoverManagers);
		DifferentialEvolution getDE(int const i);
		std::vector<InitializationType> getInitializationManagers() const;
		std::vector<MutationType> getMutationManagers() const;
		std::vector<CrossoverType> getCrossoverManagers() const;
	
		int size() const;
};