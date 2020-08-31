#pragma once
#include <vector>

#include "crossovermanager.h"
#include "solution.h"

class Particle;
class Genome : public Solution {
	public:
		Genome(int const D);
		Genome(Particle* particle);
		Genome(std::vector<double>x);
};
