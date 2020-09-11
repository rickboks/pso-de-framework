#pragma once
#include <vector>
#include "constrainthandler.h"

class Particle;

class PenaltyHandler : public ConstraintHandler {
	public:
		PenaltyHandler(std::vector<double> lb,std::vector<double> ub);
		virtual void penalize(Particle* p)=0;
};

class DeathPenalty : public PenaltyHandler {
	public:
		DeathPenalty(std::vector<double> lb,std::vector<double> ub);
		void penalize(Particle* p);
};
