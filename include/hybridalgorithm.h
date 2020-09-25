#pragma once
#include "particleupdatesettings.h"
#include "topologymanager.h"
#include "particleswarm.h"
#include "mutationmanager.h"
#include "crossovermanager.h"
//#include "selectionmanager.h"
#include "constrainthandler.h"
#include <memory>

struct HybridConfig {
	std::string const update, topology, psoCH, synchronicity, mutation, crossover, selection, adaptation, deCH;

	HybridConfig(std::string const update, std::string const topology, std::string const psoCH, std::string const synchronicity, 
			std::string const mutation, std::string const crossover, std::string const selection, 
			std::string const adaptation, std::string const deCH):
				update(update), topology(topology), psoCH(psoCH), synchronicity(synchronicity),
				mutation(mutation), crossover(crossover), selection(selection),
				adaptation(adaptation), deCH(deCH){} // For PSODE

	HybridConfig(std::string const update, std::string const topology, std::string const psoCH, std::string const synchronicity, 
			std::string const mutation, std::string const crossover, std::string const adaptation, std::string const deCH):
				update(update), topology(topology), psoCH(psoCH), synchronicity(synchronicity),
				mutation(mutation), crossover(crossover), selection(""), adaptation(adaptation), deCH(deCH){} // For PSODE2
};

class HybridAlgorithm {
	protected:
		HybridConfig const config;
	public:
		HybridAlgorithm(HybridConfig const config);
		virtual ~HybridAlgorithm() = 0;

		virtual void run(std::shared_ptr<IOHprofiler_problem<double> > const problem, 
		std::shared_ptr<IOHprofiler_csv_logger> const logger, int const evalBudget, 
		int const popSize, std::map<int,double> const particleUpdateParams) = 0;

		virtual std::string getIdString() const = 0;
};
