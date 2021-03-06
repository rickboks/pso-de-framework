#include "desuite.h"
#include "hybridsuite.h"
#include "desuite.h"
#include "psode2.h"
#include "random_suite.h"
#include <sys/types.h>

HybridAlgorithm* ha;
ParticleSwarm* pso;
DifferentialEvolution* de;
PSODE2* psode2;

void algorithm
	(std::shared_ptr<IOHprofiler_problem<double>> problem,
	std::shared_ptr<IOHprofiler_csv_logger> logger) {
	int const D = problem->IOHprofiler_get_number_of_variables(); 
	de->run(problem, logger, D*10000, 5 * D); 
	//pso->run(problem, logger, D*10000, 10 * D, {}); 
}

void _run_experiment(bool const log) {
	static registerInFactory<IOHprofiler_suite<double>,Random_suite> regSuite("random");
	de = new DifferentialEvolution(DEConfig("R1", "E", "S", "RS"));
	std::string templateFile = "./configuration.ini";

	std::string configFile = generateConfig(templateFile, de->getIdString());
	IOHprofiler_experimenter<double> experimenter(configFile,algorithm); 

	experimenter._set_independent_runs(1);
	experimenter._run();
	delete de;
}

int main(){
	bool const log = false;
	_run_experiment(log);
}
