#pragma once
#include <vector>
#include "constrainthandler.h"

class Particle;

class PenaltyHandler : public ConstraintHandler {
	public:
		PenaltyHandler(std::vector<double>const lb,std::vector<double>const ub);
		virtual void penalize(Particle* const p) const=0;
};

class DeathPenalty : public PenaltyHandler {
	public:
		DeathPenalty(std::vector<double>const lb,std::vector<double>const ub);
		void penalize(Particle* const p) const;
};
