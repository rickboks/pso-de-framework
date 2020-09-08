#include "repairhandler.h"
#include "particle.h"
#include "util.h"
#include <algorithm>
#include <limits>

RepairHandler::RepairHandler(std::vector<double>const lb, std::vector<double>const ub)
	: lb(lb),ub(ub), D(lb.size()){
}

RepairHandler::~RepairHandler(){}
void RepairHandler::repair(Particle* p){ /* do nothing */ }
void RepairHandler::repair(Particle* p, Particle* base, Particle* target){ /* do nothing */ }

// Generic
GenericRepairHandler::GenericRepairHandler(std::vector<double>const lb, std::vector<double>const ub)
	:RepairHandler(lb,ub),DEConstraintHandler(),PSOConstraintHandler(){}
GenericRepairHandler::~GenericRepairHandler(){}
void GenericRepairHandler::repair(Particle* p, Particle* base, Particle* target){}

ReinitializationRepair::ReinitializationRepair(std::vector<double>const lb, std::vector<double>const ub) : GenericRepairHandler(lb, ub){}
void ReinitializationRepair::repair(Particle* p){
	for (int i = 0; i < D; i++){
		if (p->getDimension(i) < lb[i] || p->getDimension(i) > ub[i])
			p->setDimension(i, rng.randDouble(lb[i], ub[i]));
	}
}

ProjectionRepair::ProjectionRepair(std::vector<double>const lb, std::vector<double>const ub) : GenericRepairHandler(lb, ub){}
void ProjectionRepair::repair(Particle* p){
	for (int i = 0; i < D; i++){
		if (p->getDimension(i) < lb[i]) 
			p->setDimension(i, lb[i]);
		else if (p->getDimension(i) > ub[i]) 
			p->setDimension(i, ub[i]);
	}
}

ReflectionRepair::ReflectionRepair(std::vector<double>const lb, std::vector<double>const ub) : GenericRepairHandler(lb, ub){}
void ReflectionRepair::repair(Particle* p){
	for (int i = 0; i < D; i++){
		while (true){
			if (p->getDimension(i) < lb[i]) 
				p->setDimension(i, 2 * lb[i] - p->getDimension(i));
			else if (p->getDimension(i) > ub[i]) 
				p->setDimension(i, 2 * ub[i] - p->getDimension(i));
			else 
				break;
		}
	}
}

WrappingRepair::WrappingRepair(std::vector<double>const lb, std::vector<double>const ub) : GenericRepairHandler(lb, ub){}
void WrappingRepair::repair(Particle* p){
	for (int i = 0; i < D; i++){
		while (true){
			if (p->getDimension(i) < lb[i]) 
				p->setDimension(i, ub[i] + p->getDimension(i) - lb[i]);
			else if (p->getDimension(i) > ub[i]) 
				p->setDimension(i, lb[i] + p->getDimension(i) - ub[i]);
			else 
				break;
		}
	}
}

ProjectionMidpointRepair::ProjectionMidpointRepair(std::vector<double>const lb, std::vector<double>const ub) : GenericRepairHandler(lb, ub){}
void ProjectionMidpointRepair::repair(Particle* p){
	std::vector<double> x = p->getX();
	std::vector<double>alphas(D+1);
	alphas[D] = 1;

	for (int i = 0; i < D; i++){
		if (x[i] > 0)
			alphas[i] = ub[i]/x[i];
		else if (x[i] < 0)
			alphas[i] = lb[i]/x[i];
		else
			alphas[i] = std::numeric_limits<double>::max(); //Can't divide by zero
	}

	double alpha=*std::min_element(alphas.begin(), alphas.end());
	if (alpha != 1){
		scale(x, alpha);	
		p->setX(x);
	}
}

// Differential Evolution
DERepairHandler::DERepairHandler(std::vector<double>const lb, std::vector<double>const ub): RepairHandler(lb,ub){}
DERepairHandler::~DERepairHandler(){}
void DERepairHandler::repair(Particle* p){}

RandBaseRepair::RandBaseRepair(std::vector<double>const lb, std::vector<double>const ub): DERepairHandler(lb, ub){}
void RandBaseRepair::repair(Particle* p, Particle* base, Particle* target){
	for (int i = 0; i < D; i++){
		if (p->getDimension(i) > ub[i]){
			p->setDimension(i, base->getDimension(i) + rng.randDouble(0,1) * (ub[i] - base->getDimension(i)));
		} else if (p->getDimension(i) < lb[i]){
			p->setDimension(i, base->getDimension(i) + rng.randDouble(0,1) * (lb[i] - base->getDimension(i)));
		}
	}
}

MidpointBaseRepair::MidpointBaseRepair(std::vector<double>const lb, std::vector<double>const ub): DERepairHandler(lb, ub){}
void MidpointBaseRepair::repair(Particle* p, Particle* base, Particle* target){
	for (int i = 0; i < D; i++){
		if (p->getDimension(i) > ub[i]){
			p->setDimension(i, 0.5 * (base->getDimension(i) + ub[i]));
		} else if (p->getDimension(i) < lb[i]){
			p->setDimension(i, 0.5 * (base->getDimension(i) + lb[i]));
		}
	}
}

MidpointTargetRepair::MidpointTargetRepair(std::vector<double>const lb, std::vector<double>const ub): DERepairHandler(lb, ub){}
void MidpointTargetRepair::repair(Particle* p, Particle* base, Particle* target){
	for (int i = 0; i < D; i++){
		if (p->getDimension(i) > ub[i]){
			p->setDimension(i, 0.5 * (target->getDimension(i) + ub[i]));
		} else if (p->getDimension(i) < lb[i]){
			p->setDimension(i, 0.5 * (target->getDimension(i) + lb[i]));
		}
	}
}

ConservatismRepair::ConservatismRepair(std::vector<double>const lb, std::vector<double>const ub): DERepairHandler(lb, ub){}
void ConservatismRepair::repair(Particle* p, Particle* base, Particle* target){
	for (int i = 0; i < D; i++){
		if (p->getDimension(i) > ub[i] || p->getDimension(i) < lb[i]){
			p->setDimension(i, base->getDimension(i));
		}
	}
}

//TODO resampling
