#include "penaltyhandler.h"
#include "particle.h"

PenaltyHandler::PenaltyHandler(std::vector<double> lb,std::vector<double> ub): ConstraintHandler(lb, ub){}

DeathPenalty::DeathPenalty(std::vector<double> lb,std::vector<double> ub): PenaltyHandler(lb, ub){}
void DeathPenalty::penalize(Particle* p){
	if (isFeasible(p))
		return;
	p->setFitness(std::numeric_limits<double>::max());
}
