#pragma once
#include "problem.h"
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
#include "selectionmanager.h"


class HybridAlgorithm {
	private:
		UpdateManagerType const updateManagerType;
		Topology topologyManagerType;
		std::vector<Particle*> particles;
		TopologyManager* topologyManager;
		Synchronicity const synchronicity;
		MutationType const mutationType;
		CrossoverType const crossoverType;
		SelectionType const selectionType;
		DEAdaptationType const adaptationType;
		MutationManager<Particle>* mutationManager;
		CrossoverManager<Particle>* crossoverManager;
		DEAdaptationManager* adaptationManager;
		SelectionManager* selectionManager;
		void runSynchronous(Problem const problem, int const evalBudget, int const popSize, 
			std::map<int,double> particleUpdateParams,
			double const F, double const Cr);

		void runAsynchronous(Problem const problem, int const evalBudget, int const popSize, 
			std::map<int,double> particleUpdateParams,
			double const F, double const Cr);

		public:
		HybridAlgorithm(UpdateManagerType const updateManagerType, 
			Topology topologyManager, Synchronicity const synchronous, 
			MutationType const mutationType, CrossoverType const crossoverType, 
			SelectionType const selection, DEAdaptationType adaptionType);

		~HybridAlgorithm();

		void run(Problem const problem, int const evalBudget, int const popSize, 
			std::map<int,double> particleUpdateParams, 
			double const F, double const Cr);

		void reset();
		std::string getIdString() const;

		std::vector<Particle*> copyPopulation(std::vector<Particle*>& particles);

		template <class T>
		std::vector<Solution*> toSolutions(std::vector<T*> vec);
		std::vector<Particle*> toParticles(std::vector<Solution*> sol);
	

};