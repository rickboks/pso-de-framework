#include "iohsrc/Template/Experiments/IOHprofiler_experimenter.hpp"
#include <random>
#include <set>
#include "hybridalgorithm.h"
#include "differentialevolution.h"

void genetic_algorithm
	(std::shared_ptr<IOHprofiler_problem<double>> problem,
		std::shared_ptr<IOHprofiler_csv_logger> logger) {
	int const D = problem->IOHprofiler_get_number_of_variables();
	HybridAlgorithm ha(FIPS,MULTI_SWARM, ASYNCHRONOUS, BEST_1, BINOMIAL, P3, JADE);
  DifferentialEvolution de (RANDOM, BEST_1, BINOMIAL, JADE, false);
  ParticleSwarm ps(FIPS, MULTI_SWARM, ASYNCHRONOUS);

  //ps.run(problem, logger, D*10000, D*10, std::map<int,double>());
	//ha.run(problem, logger, D*10000, D*5, std::map<int,double>(), 0.9, 0.9);
  //de.run(problem, logger, D*10000, D*10, 0.9,0.9);
}

void _run_experiment() {
  std::string configName = "./configuration.ini";
  /// An example for PBO suite.
  IOHprofiler_experimenter<double> experimenter(configName,"PSO",genetic_algorithm); 
  /// An exmaple for BBOB suite.
  /// IOHprofiler_experimenter<double> experimenter(configName, random_search);
  experimenter._set_independent_runs(10);
  experimenter._run();
}

int main(){
  _run_experiment();
  return 0;
}