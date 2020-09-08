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

void RepairHandler::repairVelocity(Particle* p, int i){
	if (p->isPSO){ // Only repair is we're dealing with a 'PSO particle'
		p->setV(i, -0.5 * p->getV(i)); // Deterministic back velocity repair from paper pso2.pdf
	}
}

// Generic
GenericRepairHandler::GenericRepairHandler(std::vector<double>const lb, std::vector<double>const ub)
	:RepairHandler(lb,ub),DEConstraintHandler(),PSOConstraintHandler(){}
GenericRepairHandler::~GenericRepairHandler(){}
void GenericRepairHandler::repair(Particle* p, Particle* base, Particle* target){}

ReinitializationRepair::ReinitializationRepair(std::vector<double>const lb, std::vector<double>const ub) : GenericRepairHandler(lb, ub){}
void ReinitializationRepair::repair(Particle* p){
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i] || p->getX(i) > ub[i]){
			p->setX(i, rng.randDouble(lb[i], ub[i]));
			repairVelocity(p, i);
		}
	}
}

ProjectionRepair::ProjectionRepair(std::vector<double>const lb, std::vector<double>const ub) : GenericRepairHandler(lb, ub){}
void ProjectionRepair::repair(Particle* p){
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i]){
			p->setX(i, lb[i]);
			repairVelocity(p, i);
		} else if (p->getX(i) > ub[i]){
			p->setX(i, ub[i]);
			repairVelocity(p, i);
		}
	}
}

ReflectionRepair::ReflectionRepair(std::vector<double>const lb, std::vector<double>const ub) : GenericRepairHandler(lb, ub){}
void ReflectionRepair::repair(Particle* p){
	for (int i = 0; i < D; i++){
		bool is_repaired = false;
		while (p->getX(i) < lb[i]){
			p->setX(i, 2 * lb[i] - p->getX(i));
			is_repaired = true;
		}
		while (p->getX(i) > ub[i]){
			p->setX(i, 2 * ub[i] - p->getX(i));
			is_repaired = true;
		}
		if (is_repaired)
			repairVelocity(p, i);
	}
}

WrappingRepair::WrappingRepair(std::vector<double>const lb, std::vector<double>const ub) : GenericRepairHandler(lb, ub){}
void WrappingRepair::repair(Particle* p){
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i]){
			p->setX(i, ub[i] - std::fmod(lb[i] - p->getX(i), std::fabs(ub[i]-lb[i])));
			repairVelocity(p,i);
		} else if (p->getX(i) > ub[i]) {
			p->setX(i, lb[i] - std::fmod(p->getX(i) - ub[i], std::fabs(ub[i]-lb[i])));
			repairVelocity(p,i);
		}
	}
}

ProjectionMidpointRepair::ProjectionMidpointRepair(std::vector<double>const lb, std::vector<double>const ub) : GenericRepairHandler(lb, ub){}
void ProjectionMidpointRepair::repair(Particle* p){
	std::vector<double> x = p->getX();
	std::vector<double>alphas(D+1);
	alphas[D] = 1;

	for (int i = 0; i < D; i++){
		if (x[i] > 0){
			alphas[i] = ub[i]/x[i];
			repairVelocity(p,i);
		} else if (x[i] < 0){
			alphas[i] = lb[i]/x[i];
			repairVelocity(p,i);
		} else
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
		if (p->getX(i) > ub[i]){
			p->setX(i, base->getX(i) + rng.randDouble(0,1) * (ub[i] - base->getX(i)));
		} else if (p->getX(i) < lb[i]){
			p->setX(i, base->getX(i) + rng.randDouble(0,1) * (lb[i] - base->getX(i)));
		}
	}
}

MidpointBaseRepair::MidpointBaseRepair(std::vector<double>const lb, std::vector<double>const ub): DERepairHandler(lb, ub){}
void MidpointBaseRepair::repair(Particle* p, Particle* base, Particle* target){
	for (int i = 0; i < D; i++){
		if (p->getX(i) > ub[i]){
			p->setX(i, 0.5 * (base->getX(i) + ub[i]));
		} else if (p->getX(i) < lb[i]){
			p->setX(i, 0.5 * (base->getX(i) + lb[i]));
		}
	}
}

MidpointTargetRepair::MidpointTargetRepair(std::vector<double>const lb, std::vector<double>const ub): DERepairHandler(lb, ub){}
void MidpointTargetRepair::repair(Particle* p, Particle* base, Particle* target){
	for (int i = 0; i < D; i++){
		if (p->getX(i) > ub[i]){
			p->setX(i, 0.5 * (target->getX(i) + ub[i]));
		} else if (p->getX(i) < lb[i]){
			p->setX(i, 0.5 * (target->getX(i) + lb[i]));
		}
	}
}

ConservatismRepair::ConservatismRepair(std::vector<double>const lb, std::vector<double>const ub): DERepairHandler(lb, ub){}
void ConservatismRepair::repair(Particle* p, Particle* base, Particle* target){
	for (int i = 0; i < D; i++){
		if (p->getX(i) > ub[i] || p->getX(i) < lb[i]){
			p->setX(i, base->getX(i));
		}
	}
}
 
//TODO resampling
