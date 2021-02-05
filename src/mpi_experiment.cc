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
#include "random_suite.h"

DESuite suite;

void experiment
	(std::shared_ptr<IOHprofiler_problem<double>> problem,
		std::shared_ptr<IOHprofiler_csv_logger> logger) {

	int const D = problem->IOHprofiler_get_number_of_variables();
	int const popSize = 100;

	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	DifferentialEvolution de = suite.getDE(id);
  	de.run(problem, logger, D*10000, popSize);
}

int main(int argc, char **argv) {
	static registerInFactory<IOHprofiler_suite<double>,Random_suite> regSuite("random");
	suite.setDEAdaptationManagers({"S"});
	suite.setCrossoverManagers({"B"});

	std::string const templateFile = "./configuration.ini";
	int id;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	if (id < suite.size()){
		std::string const configFile = generateConfig(templateFile, suite.getDE(id).getIdString());
		IOHprofiler_experimenter<double> experimenter(configFile, experiment);
		experimenter._set_independent_runs(100);
		experimenter._run();
	} else {
		std::cerr << "Error: suite does not contain " << id << std::endl;
	}

	MPI_Finalize();

	return 0;
}
