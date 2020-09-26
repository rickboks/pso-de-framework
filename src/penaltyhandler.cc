#include "penaltyhandler.h"
#include "particle.h"

DeathPenalty::DeathPenalty(std::vector<double>const lb,std::vector<double>const ub)
	: ConstraintHandler(lb,ub), DEConstraintHandler(lb,ub), PSOConstraintHandler(lb,ub){}

void DeathPenalty::penalize(Solution* const p) const {
	if (!isFeasible(p))
		p->setFitness(std::numeric_limits<double>::max());
}
