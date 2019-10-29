#include "iohsrc/Template/Experiments/IOHprofiler_experimenter.hpp"
#include <random>
#include <set>
#include "hybridalgorithm.h"
#include "differentialevolution.h"

HybridAlgorithm* ha;
ParticleSwarm* pso;
DifferentialEvolution* de;

std::vector<MutationType> mutations ({
  BEST_1,
  BEST_2,
  RAND_1,
  RAND_2,
  TTB_1,
  TTPB_1, 
});

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
      //for (int i = 0; i < MUT_END; i++){
      ha = new HybridAlgorithm(DECR_INERTIA_WEIGHT, VON_NEUMANN, ASYNCHRONOUS, TTPB_1, BINOMIAL, P3, JADE);
      //std::cout << ha->getIdString() << std::endl;
      //de = new DifferentialEvolution(OPPOSITION, BEST_1, BINOMIAL, JADE, false);

      std::string configName = "./configuration.ini";
      IOHprofiler_experimenter<double> experimenter(configName,pso->getIdString(),algorithm); 
      experimenter._set_independent_runs(10);
      experimenter._run();
      delete ha;
      delete pso;
}
  //}


int main(){
  _run_experiment();
  return 0; 
}