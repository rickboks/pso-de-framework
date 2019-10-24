#pragma once
#include "problem.h"
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
		Problem problem;

		void initRandom(std::vector<Genome*>& genomes) const;
		void initOpposition(std::vector<Genome*>& genomes) const;

	public:
		DEInitializer(DEInitializationType const type, Problem problem);
		void initialize(std::vector<Genome*>& genomes) const;
};