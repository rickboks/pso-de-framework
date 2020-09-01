#include <IOHprofiler_experimenter.h>
#include <random>
#include <set>
#include "hybridalgorithm.h"
#include "differentialevolution.h"
#include "hybridsuite.h"
#include "psode2.h"
#include "psode.h"
HybridAlgorithm* ha;
ParticleSwarm* pso;
DifferentialEvolution* de;
PSODE* psode;
PSODE2* psode2;

void algorithm
(std::shared_ptr<IOHprofiler_problem<double>> problem,
 std::shared_ptr<IOHprofiler_csv_logger> logger) {
    int const D = problem->IOHprofiler_get_number_of_variables(); 
    //psode2->run(problem, logger, D*10000, D*5, std::map<int,double>()); 
    //psode->run(problem, logger, D*10000, D*5, std::map<int,double>()); 
    de->run(problem, logger, D*10000, D*5); 
    //pso->run(problem, logger, D*10000, D*5, std::map<int,double>()); 
}

void _run_experiment(bool log) {
    psode2 = new PSODE2(DECR_INERTIA_WEIGHT, VON_NEUMANN, ASYNCHRONOUS, TTB_1, BINOMIAL, P3, JADE);
    psode = new PSODE(INERTIA_WEIGHT, VON_NEUMANN, ASYNCHRONOUS, TTB_1, BINOMIAL, P3, JADE);
    pso = new ParticleSwarm(INERTIA_WEIGHT, VON_NEUMANN, ASYNCHRONOUS);
	de = new DifferentialEvolution(RANDOM, TTB_1, BINOMIAL, JADE, false);

    if (log)
		pso->enableLogging();

    std::string configName = "./configuration.ini";
    IOHprofiler_experimenter<double> experimenter(configName,algorithm); 
    experimenter._set_independent_runs(5);
    experimenter._run();
    delete pso;
	delete de;
	delete psode;
	delete psode2;
}

int main(int argc, char** argv){
    bool log = false;
    _run_experiment(log);	
}
