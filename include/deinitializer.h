#pragma once
#include "problem.h"
#include <vector>
class Genome;

typedef void (*evaluate_function_t)(const double *x, double *y);
enum InitializationType {
	RANDOM,
	OPPOSITION,
	INIT_END
};

class DEInitializer {
	private:
		InitializationType const initializationType;
		Problem problem;

		void initRandom(std::vector<Genome*>& genomes) const;
		void initOpposition(std::vector<Genome*>& genomes) const;

	public:
		DEInitializer(InitializationType const type, Problem problem);
		void initialize(std::vector<Genome*>& genomes) const;
};