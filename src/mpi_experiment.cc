#include "iohsrc/Template/Experiments/IOHprofiler_experimenter.hpp"
#include <random>
#include <set>
#include <mpi.h>
#include "hybridalgorithm.h"
#include "differentialevolution.h"
#include "desuite.h"
#include "hybridsuite.h"
#include "particleswarmsuite.h"

std::vector<UpdateManagerType> updateManagers ({
	FIPS,
	BARE_BONES,
	INERTIA_WEIGHT,
	DECR_INERTIA_WEIGHT
});

std::vector<Topology> topologies ({
	LBEST,
	GBEST,
	VON_NEUMANN,
	MULTI_SWARM,	
	INCREASING
});

std::vector<MutationType> mutations ({
	BEST_1,
	BEST_2,
	TTB_1,
	TTPB_1,
	TO1
});

std::vector<CrossoverType> crossovers ({
	BINOMIAL,
	EXPONENTIAL
});

std::vector<DEAdaptationType> adaptationTypes ({
	JADE
});

std::vector<SelectionType> selectionTypes({
	P2, P3, U2, U3
});

std::vector<Synchronicity> synchronicities({
	ASYNCHRONOUS
});

HybridSuite suite;

void experiment
	(std::shared_ptr<IOHprofiler_problem<double>> problem,
		std::shared_ptr<IOHprofiler_csv_logger> logger) {
	int const D = problem->IOHprofiler_get_number_of_variables();
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	HybridAlgorithm h = suite.getHybrid(id);
  	h.run(problem, logger, D*10000, D*5, std::map<int,double>());
}

int main(int argc, char **argv) {
	suite.setUpdateManagers(updateManagers);
	suite.setTopologyManagers(topologies);
	suite.setMutationManagers(mutations);
	suite.setCrossoverManagers(crossovers);
	suite.setSelectionManagers(selectionTypes);
	suite.setSynchronicities(synchronicities);
	suite.setDEAdaptationManagers(adaptationTypes);

	std::string configFile = "./configuration.ini";
	int id;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	if (id < suite.size()){
		IOHprofiler_experimenter<double> experimenter(configFile,suite.getHybrid(id).getIdString(), experiment , true);
	  	experimenter._set_independent_runs(30);
	  	experimenter._run();
	} else {
		std::cerr << "Error: suite does not contain " << id << std::endl;
	}
	MPI_Finalize();

	return 0;
}
