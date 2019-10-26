#include "iohsrc/Template/Experiments/IOHprofiler_experimenter.hpp"
#include <random>
#include <set>
#include "hybridalgorithm.h"
#include "differentialevolution.h"

void algorithm
	(std::shared_ptr<IOHprofiler_problem<double>> problem,
		std::shared_ptr<IOHprofiler_csv_logger> logger) {
	int const D = problem->IOHprofiler_get_number_of_variables();
  HybridAlgorithm ha(DECR_INERTIA_WEIGHT, VON_NEUMANN, ASYNCHRONOUS, BEST_1, BINOMIAL, U2, JADE);
	ha.run(problem, logger, D*10000, D*5, std::map<int,double>(), 0.9, 0.9);  
}

void _run_experiment() {
    std::string configName = "./configuration.ini";
    IOHprofiler_experimenter<double> experimenter(configName,"Hybrid",algorithm); 
    experimenter._set_independent_runs(10);
    experimenter._run();
  }


int main(){
  _run_experiment();
  return 0;
}