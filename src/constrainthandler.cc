#include "constrainthandler.h"
#include "particle.h"
#include "repairhandler.h"
#include "penaltyhandler.h"
#define LC(X) [](std::vector<double>lb, std::vector<double>ub){return new X(lb,ub);}

ConstraintHandler::ConstraintHandler(std::vector<double>const lb,std::vector<double>const ub)
: lb(lb), ub(ub), D(lb.size()){}

ConstraintHandler::~ConstraintHandler(){}

void ConstraintHandler::repair(Particle* const p, Particle const* const base, Particle const* const target) const{}

void ConstraintHandler::repair(Particle* const p) const{}

bool ConstraintHandler::resample(Particle const * const p, int const resamples) const{
	return false;
}
void ConstraintHandler::repairVelocityPre(Particle* const p) const{}

void ConstraintHandler::repairVelocityPost(Particle* const p, int const i) const{
	if (p->isPSO) p->setV(i, -0.5 * p->getV(i));
}

bool ConstraintHandler::isFeasible(Particle const * const p) const{
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i] || p->getX(i) > ub[i]){
			return false;
		}
	}
	return true;
}

std::map<std::string, std::function<ConstraintHandler* (std::vector<double>, std::vector<double>)>> const psoCHs ({
		// Generic
		{"reinitialization", LC(ReinitializationRepair)},
		{"projection", LC(ProjectionRepair)},
		{"reflection", LC(ReflectionRepair)},
		{"wrapping", LC(WrappingRepair)},
		{"death penalty", LC(DeathPenalty)},
		{"resampling", LC(ResamplingRepair)},

		// PSO
		{"hyperbolic", LC(HyperbolicRepair)},
		{"pbest", LC(PBestDimRepair)},
});

std::map<std::string, std::function<ConstraintHandler* (std::vector<double>, std::vector<double>)>> const deCHs ({
		// Generic
		{"reinitialization", LC(ReinitializationRepair)},
		{"projection", LC(ProjectionRepair)},
		{"reflection", LC(ReflectionRepair)},
		{"wrapping", LC(WrappingRepair)},
		{"death penalty", LC(DeathPenalty)},
		{"resampling", LC(ResamplingRepair)},

		// DE
		{"rand base", LC(RandBaseRepair)},
		{"midpoint base", LC(MidpointBaseRepair)},
		{"midpoint target", LC(MidpointTargetRepair)},
		{"conservatism", LC(ConservatismRepair)},
		{"projection midpoint", LC(ProjectionMidpointRepair)},
		{"projection base", LC(ProjectionBaseRepair)},
});

//GenericConstraintHandler::GenericConstraintHandler(std::vector<double> lb,std::vector<double> ub)
//: ConstraintHandler(lb,ub){}

//DEConstraintHandler::DEConstraintHandler(std::vector<double> lb,std::vector<double> ub)
//: ConstraintHandler(lb,ub){}

//PSOConstraintHandler::PSOConstraintHandler(std::vector<double> lb,std::vector<double> ub)
//: ConstraintHandler(lb,ub){}


//DEConstraintHandler::~DEConstraintHandler(){}

//PSOConstraintHandler::~PSOConstraintHandler(){}
