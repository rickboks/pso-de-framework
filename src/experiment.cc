#include "iohsrc/Template/Experiments/IOHprofiler_experimenter.hpp"
#include <random>
#include <set>
#include "hybridalgorithm.h"
#include "differentialevolution.h"
HybridAlgorithm* ha;
ParticleSwarm* pso;
DifferentialEvolution* de;

//5 dimensions for HA

void algorithm
	(std::shared_ptr<IOHprofiler_problem<double>> problem,
		std::shared_ptr<IOHprofiler_csv_logger> logger) {
	int const D = problem->IOHprofiler_get_number_of_variables();

	//ha->run(problem, logger, D*10000, D*5, std::map<int,double>());  
  //de->run(problem, logger, D*10000, D*5);
  pso->run(problem, logger, D*10000, D*5, std::map<int,double>());
}


void _run_experiment() {
    pso = new ParticleSwarm(DECR_INERTIA_WEIGHT, VON_NEUMANN, ASYNCHRONOUS);
    ha = new HybridAlgorithm(DECR_INERTIA_WEIGHT, VON_NEUMANN, ASYNCHRONOUS, BEST_1, BINOMIAL, P3, JADE);
    de = new DifferentialEvolution(OPPOSITION, BEST_1, BINOMIAL, JADE, false);

    std::string configName = "./configuration.ini";
    IOHprofiler_experimenter<double> experimenter(configName,pso->getIdString(),algorithm); 
    experimenter._set_independent_runs(10);
    experimenter._run();

    delete pso;
    delete ha;
    delete de;
  }


int main(){
  _run_experiment();
  return 0;
}