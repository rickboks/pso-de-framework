#include "penaltyhandler.h"
#include "particle.h"

PenaltyHandler::PenaltyHandler(std::vector<double>const lb,std::vector<double>const ub): ConstraintHandler(lb, ub){}

DeathPenalty::DeathPenalty(std::vector<double>const lb,std::vector<double>const ub): PenaltyHandler(lb, ub){}

void DeathPenalty::penalize(Particle* const p) const {
	if (!isFeasible(p))
		p->setFitness(std::numeric_limits<double>::max());
}
