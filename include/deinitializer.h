#pragma once
#include "solution.h"
#include <vector>
class Genome;

typedef void (*evaluate_function_t)(const double *x, double *y);
enum DEInitializationType {
	RANDOM,
	OPPOSITION,
	INIT_END
};

class DEInitializer {
	private:
		DEInitializationType const initializationType;
		std::shared_ptr<IOHprofiler_problem<double> > problem;
		std::shared_ptr<IOHprofiler_csv_logger> logger;

		void initRandom(std::vector<Genome*>& genomes) const;
		void initOpposition(std::vector<Genome*>& genomes) const;

	public:
		DEInitializer(DEInitializationType const type, std::shared_ptr<IOHprofiler_problem<double> > problem,std::shared_ptr<IOHprofiler_csv_logger> logger);
		void initialize(std::vector<Genome*>& genomes) const;
};