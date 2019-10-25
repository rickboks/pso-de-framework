#include "iohsrc/Template/Experiments/IOHprofiler_experimenter.hpp"
#include <random>
#include <set>
#include <mpi.h>
#include "hybridalgorithm.h"
#include "differentialevolution.h"
#include "desuite.h"
#include "hybridsuite.h"

HybridSuite s;
void experiment
	(std::shared_ptr<IOHprofiler_problem<double>> problem,
		std::shared_ptr<IOHprofiler_csv_logger> logger) {
	int const D = problem->IOHprofiler_get_number_of_variables();
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);



	HybridAlgorithm h = s.getHybrid(id);
  	h.run(problem, logger, D*10000, D*10, std::map<int,double>(), 0.9,0.9);
}

int main(int argc, char **argv) {
	std::string configName = "./configuration.ini";
	int id;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	if (id < s.size()){
		IOHprofiler_experimenter<double> experimenter(configName,s.getHybrid(id).getIdString(), experiment);
	  	experimenter._set_independent_runs(10);
	  	experimenter._run();
	} else {
		std::cout << "Error: suite does not contain " << id << std::endl;
	}
	MPI_Finalize();

	return 0;
}



void _run_experiment() {

}