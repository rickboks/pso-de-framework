#include "constrainthandler.h"
#include "particle.h"
#include "repairhandler.h"
#include "penaltyhandler.h"
#define LC(X) [](std::vector<double>lb, std::vector<double>ub){return new X(lb,ub);}

bool ConstraintHandler::isFeasible(Solution const * const p) const{
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i] || p->getX(i) > ub[i]){
			return false;
		}
	}
	return true;
}

bool ConstraintHandler::resample(Solution const * const p, int const resamples){
	return false;
}

std::map<std::string, std::function<DEConstraintHandler*(std::vector<double>, std::vector<double>)>> const deCHs ({
	// Generic
	{"DP", LC(DeathPenalty)},
	{"RS", LC(ResamplingRepair)},

	// Almost generic
	{"RI", LC(DEReinitializationRepair)},
	{"PR", LC(DEProjectionRepair)},
	{"RF", LC(DEReflectionRepair)},
	{"WR", LC(DEWrappingRepair)},
	{"TR", LC(DETransformationRepair)},

	//// DE
	{"RB", LC(RandBaseRepair)},
	{"MB", LC(MidpointBaseRepair)},
	{"MT", LC(MidpointTargetRepair)},
	{"CO", LC(ConservatismRepair)},
	{"PM", LC(ProjectionMidpointRepair)},
	{"PB", LC(ProjectionBaseRepair)},
});

std::map<std::string, std::function<PSOConstraintHandler*(std::vector<double>, std::vector<double>)>> const psoCHs {
	// Generic
	{"DP", LC(DeathPenalty)},
	{"RS", LC(ResamplingRepair)},

	// Almost generic
	{"RI", LC(PSOReinitializationRepair)},
	{"PR", LC(PSOProjectionRepair)},
	{"RF", LC(PSOReflectionRepair)},
	{"WR", LC(PSOWrappingRepair)},
	{"TR", LC(PSOTransformationRepair)},

	//// PSO
	{"HY", LC(HyperbolicRepair)},
	{"PD", LC(PBestDimRepair)},
};
