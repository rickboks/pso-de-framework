#pragma once
#include <vector>

#include "crossovermanager.h"
#include "solution.h"

class Particle;
class Genome : public Solution {
	private:
		double weightFactor; //used only for population-based mutation
	public:
		Genome(int const D);
		Genome(Particle* particle);
		Genome(std::vector<double>x);
		void setWeightFactor(double weightFactor);
		double getWeightFactor() const;
};