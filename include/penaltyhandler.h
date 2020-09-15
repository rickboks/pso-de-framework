#pragma once
#include <vector>
#include "constrainthandler.h"

class Particle;

class PenaltyHandler : public ConstraintHandler {
	public:
		PenaltyHandler(std::vector<double>const lb,std::vector<double>const ub);
};

class DeathPenalty : public PenaltyHandler {
	public:
		DeathPenalty(std::vector<double>const lb,std::vector<double>const ub);
		void penalize(Particle* const p) const;
};
