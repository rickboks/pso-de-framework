#include "constrainthandler.h"
#include "particle.h"
#include "repairhandler.h"
#define LC(X) [](std::vector<double>lb, std::vector<double>ub){return new X(lb,ub);}

bool ConstraintHandler::isFeasible(Solution const * const p) const{
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i] - 1.0e-12 || p->getX(i) > ub[i] + 1.0e-12){
			return false;
		}
	}
	return true;
}

bool ConstraintHandler::resample(Solution * const p, int const resamples){
	return false;
}

int ConstraintHandler::getCorrections() const {
	return nCorrected;
}

std::map<std::string, std::function<DEConstraintHandler*(std::vector<double>, std::vector<double>)>> const deCHs ({
	// Generic
	{"DP", LC(DeathPenalty)},
	{"RS", LC(ResamplingRepair)},

	// Almost generic
	{"RI", LC(ReinitializationRepair)},
	{"PR", LC(ProjectionRepair)},
	{"RF", LC(ReflectionRepair)},
	{"WR", LC(WrappingRepair)},
	{"TR", LC(TransformationRepair)},

	//// DE
	{"RB", LC(RandBaseRepair)},
	{"MB", LC(MidpointBaseRepair)},
	{"MT", LC(MidpointTargetRepair)},
	{"PM", LC(ProjectionMidpointRepair)},
	{"PB", LC(ProjectionBaseRepair)},
	{"CO", LC(ConservatismRepair)}
});

std::map<std::string, std::function<PSOConstraintHandler*(std::vector<double>, std::vector<double>)>> const psoCHs {
	// Generic
	{"DP", LC(DeathPenalty)},
	{"RS", LC(ResamplingRepair)},

	// Almost generic
	{"RI", LC(ReinitializationRepair)},
	{"PR", LC(ProjectionRepair)},
	{"RF", LC(ReflectionRepair)},
	{"WR", LC(WrappingRepair)},
	{"TR", LC(TransformationRepair)},

	//// PSO
	{"HY", LC(HyperbolicRepair)},
	{"PD", LC(PBestDimRepair)},
};
