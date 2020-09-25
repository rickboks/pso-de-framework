#include "constrainthandler.h"
#include "particle.h"
#include "repairhandler.h"
#include "penaltyhandler.h"
#define LC(X) [](std::vector<double>lb, std::vector<double>ub){return new X(lb,ub);}

ConstraintHandler::ConstraintHandler(std::vector<double>const lb,std::vector<double>const ub): lb(lb), ub(ub), D(lb.size()){
}
ConstraintHandler::~ConstraintHandler(){};

bool ConstraintHandler::isFeasible(Solution const * const p) const{
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i] || p->getX(i) > ub[i]){
			return false;
		}
	}
	return true;
}

bool ConstraintHandler::resample(Solution const * const p, int const resamples) const{
	return false;
}

DEConstraintHandler::DEConstraintHandler(std::vector<double>const lb,std::vector<double>const ub)
: ConstraintHandler(lb, ub){}

DEConstraintHandler::~DEConstraintHandler(){}


PSOConstraintHandler::PSOConstraintHandler(std::vector<double>const lb,std::vector<double>const ub)
: ConstraintHandler(lb, ub){}
PSOConstraintHandler::~PSOConstraintHandler(){};
//void DEConstraintHandler::repairVelocityPre(Solution* const p) const{}
//void DEConstraintHandler::repairVelocityPost(Solution* const p, int const i) const{
	////if (p->isPSO) 
		////p->setV(i, -0.5 * p->getV(i));
//}

//std::map<std::string, std::function<ConstraintHandler* (std::vector<double>, std::vector<double>)>> const psoCHs ({
	//// Generic
	//{"RI", LC(ReinitializationRepair)},
	//{"PR", LC(ProjectionRepair)},
	//{"RF", LC(ReflectionRepair)},
	//{"WR", LC(WrappingRepair)},
	//{"DP", LC(DeathPenalty)},
	//{"RS", LC(ResamplingRepair)},
	//{"TR", LC(TransformationRepair)},

	//// PSO
	//{"HY", LC(HyperbolicRepair)},
	//{"PD", LC(PBestDimRepair)},
//});

std::map<std::string, std::function<ConstraintHandler*(std::vector<double>, std::vector<double>)>> const deCHs ({
	// Generic
	{"RI", LC(ReinitializationRepair)},
	{"PR", LC(ProjectionRepair)},
	{"RF", LC(ReflectionRepair)},
	{"WR", LC(WrappingRepair)},
	{"DP", LC(DeathPenalty)},
	{"RS", LC(ResamplingRepair)},
	{"TR", LC(TransformationRepair)},

	//// DE
	{"RB", LC(RandBaseRepair)},
	{"MB", LC(MidpointBaseRepair)},
	{"MT", LC(MidpointTargetRepair)},
	{"CO", LC(ConservatismRepair)},
	{"PM", LC(ProjectionMidpointRepair)},
	{"PB", LC(ProjectionBaseRepair)},
});

std::map<std::string, std::function<ConstraintHandler*(std::vector<double>, std::vector<double>)>> const psoCHs {
	// Generic
	{"RI", LC(ReinitializationRepair)},
	{"PR", LC(ProjectionRepair)},
	{"RF", LC(ReflectionRepair)},
	{"WR", LC(WrappingRepair)},
	{"DP", LC(DeathPenalty)},
	{"RS", LC(ResamplingRepair)},
	{"TR", LC(TransformationRepair)},

	//// PSO
	{"HY", LC(HyperbolicRepair)},
	{"PD", LC(PBestDimRepair)},
};

//GenericConstraintHandler::GenericConstraintHandler(std::vector<double> lb,std::vector<double> ub)
//: ConstraintHandler(lb,ub){}

//DEConstraintHandler::DEConstraintHandler(std::vector<double> lb,std::vector<double> ub)
//: ConstraintHandler(lb,ub){}

//PSOConstraintHandler::PSOConstraintHandler(std::vector<double> lb,std::vector<double> ub)
//: ConstraintHandler(lb,ub){}

//DEConstraintHandler::~DEConstraintHandler(){}

//PSOConstraintHandler::~PSOConstraintHandler(){}
