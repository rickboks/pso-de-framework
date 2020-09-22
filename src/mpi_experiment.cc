#include <IOHprofiler_experimenter.h>
#include <random>
#include <set>
#include <mpi.h>
#include <fstream>
#include "hybridalgorithm.h"
#include "differentialevolution.h"
#include "desuite.h"
#include "hybridsuite.h"
#include "particleswarmsuite.h"
#include "util.h"

DESuite suite;

void experiment
	(std::shared_ptr<IOHprofiler_problem<double>> problem,
		std::shared_ptr<IOHprofiler_csv_logger> logger) {
	int const D = problem->IOHprofiler_get_number_of_variables();
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	DifferentialEvolution de = suite.getDE(id);
  	de.run(problem, logger, D*10000, D*5);
}

int main(int argc, char **argv) {
	suite.setMutationManagers(std::vector<std::string>{
		"R1", "R2", "B1", "T1"
	});

	suite.setCrossoverManagers(std::vector<std::string>{
		"B", "E"
	});

	suite.setDEAdaptationManagers(std::vector<std::string>{
		"N"
	});

	suite.setConstraintHandlers(std::vector<std::string>{
		"PR", "WR", "DP"
	});

	std::string templateFile = "./configuration.ini";
	int id;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	if (id < suite.size()){
		std::string configFile = generateConfig(templateFile, suite.getDE(id).getIdString());
		std::cout << configFile << std::endl;
		IOHprofiler_experimenter<double> experimenter(configFile, experiment);
		experimenter._set_independent_runs(10);
		experimenter._run();
	} else {
		std::cerr << "Error: suite does not contain " << id << std::endl;
	}

	MPI_Finalize();

	return 0;
}
