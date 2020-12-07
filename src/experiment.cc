#include <IOHprofiler_experimenter.h>
#include <random>
#include <set>
#include "desuite.h"
#include "hybridalgorithm.h"
#include "differentialevolution.h"
#include "hybridsuite.h"
#include "psode2.h"
#include "util.h"

HybridAlgorithm* ha;
ParticleSwarm* pso;
DifferentialEvolution* de;
PSODE2* psode2;

void algorithm
(std::shared_ptr<IOHprofiler_problem<double>> problem,
 std::shared_ptr<IOHprofiler_csv_logger> logger) {
    int const D = problem->IOHprofiler_get_number_of_variables(); 
    //psode2->run(problem, logger, D*10000, D*5, std::map<int,double>()); 
    //psode->run(problem, logger, D*10000, D*5, std::map<int,double>()); 
    //pso->run(problem, logger, D*10000, D*5, std::map<int,double>()); 
    de->run(problem, logger, D*10000, 5 * D); 
}

void _run_experiment(bool const log) {

    //DESuite suite;
    //std::cout << "size: " << suite.size() << std::endl;
    //return;
    //psode2 = new PSODE2(HybridConfig("I", "N", "HY", "A", "T1", "B", "J", "RI"));
    //psode = new PSODE(HybridConfig("I", "N", "HY", "A", "T1", "B", "P3", "J", "RI"));    
    //pso = new ParticleSwarm(PSOConfig("I", "M", "PR", "A"));

    de = new DifferentialEvolution(DEConfig("P1", "B", "S", "PM"));
	std::string templateFile = "./configuration.ini";
    std::string configFile = generateConfig(templateFile, de->getIdString());
    IOHprofiler_experimenter<double> experimenter(configFile,algorithm); 

    experimenter._set_independent_runs(5);
    experimenter._run();
    delete de;
    //delete de;
	//delete psode;
	//delete psode2;
}

int main(){
    bool const log = false;
    _run_experiment(log);
}
