#include "repairhandler.h"
#include "particle.h"
#include "util.h"
#include <algorithm>
#include <limits>

RepairHandler::~RepairHandler(){}
RepairHandler::RepairHandler(std::vector<double>const lb, std::vector<double>const ub) : ConstraintHandler(lb, ub){}
void RepairHandler::repair(Particle* const p){}
void RepairHandler::repair(Particle* const p, Particle* const base, Particle* const target){}

ReinitializationRepair::ReinitializationRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void ReinitializationRepair::repair(Particle* const p){
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i] || p->getX(i) > ub[i]){
			p->setX(i, rng.randDouble(lb[i], ub[i]));
			repairVelocity(p, i);
		}
	}
}

ProjectionRepair::ProjectionRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void ProjectionRepair::repair(Particle* const p){
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

ReflectionRepair::ReflectionRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void ReflectionRepair::repair(Particle* const p){
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

WrappingRepair::WrappingRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void WrappingRepair::repair(Particle* const p){
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

// Differential Evolution
RandBaseRepair::RandBaseRepair(std::vector<double>const lb, std::vector<double>const ub): RepairHandler(lb, ub){}
void RandBaseRepair::repair(Particle* const p, Particle* const base, Particle* const target){
	for (int i = 0; i < D; i++){
		if (p->getX(i) > ub[i]){
			p->setX(i, base->getX(i) + rng.randDouble(0,1) * (ub[i] - base->getX(i)));
		} else if (p->getX(i) < lb[i]){
			p->setX(i, base->getX(i) + rng.randDouble(0,1) * (lb[i] - base->getX(i)));
		}
	}
}

MidpointBaseRepair::MidpointBaseRepair(std::vector<double>const lb, std::vector<double>const ub): RepairHandler(lb, ub){}
void MidpointBaseRepair::repair(Particle* const p, Particle* const base, Particle* const target){
	for (int i = 0; i < D; i++){
		if (p->getX(i) > ub[i]){
			p->setX(i, 0.5 * (base->getX(i) + ub[i]));
		} else if (p->getX(i) < lb[i]){
			p->setX(i, 0.5 * (base->getX(i) + lb[i]));
		}
	}
}

MidpointTargetRepair::MidpointTargetRepair(std::vector<double>const lb, std::vector<double>const ub): RepairHandler(lb, ub){}
void MidpointTargetRepair::repair(Particle* const p, Particle* const base, Particle* const target){
	for (int i = 0; i < D; i++){
		if (p->getX(i) > ub[i]){
			p->setX(i, 0.5 * (target->getX(i) + ub[i]));
		} else if (p->getX(i) < lb[i]){
			p->setX(i, 0.5 * (target->getX(i) + lb[i]));
		}
	}
}

ConservatismRepair::ConservatismRepair(std::vector<double>const lb, std::vector<double>const ub): RepairHandler(lb, ub){}
void ConservatismRepair::repair(Particle* const p, Particle* const base, Particle* const target){
	for (int i = 0; i < D; i++){
		if (p->getX(i) > ub[i] || p->getX(i) < lb[i]){
			p->setX(i, base->getX(i));
		}
	}
}

ProjectionMidpointRepair::ProjectionMidpointRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void ProjectionMidpointRepair::repair(Particle* const p, Particle* const base, Particle* const target){
	std::vector<double> x = p->getX();
	std::vector<double>alphas(D+1);
	alphas[D] = 1;

	for (int i = 0; i < D; i++){
		if (x[i] > ub[i]){
			alphas[i] = (lb[i] - ub[i])/(lb[i] - 2 * x[i] + ub[i]);
		} else if (x[i] < lb[i]){
			alphas[i] = (ub[i] - lb[i])/(lb[i] - 2 * x[i] + ub[i]);
		} else
			alphas[i] = std::numeric_limits<double>::max(); //Can't divide by zero
	}

	double alpha=*std::min_element(alphas.begin(), alphas.end());
	if (alpha != 1){
		std::vector<double> midpoint(D);
		add(lb, ub, midpoint);
		scale(midpoint, 0.5*(1-alpha));
		scale(x, alpha);
		add(x, midpoint, x);
		p->setX(x);
	}

}

ResamplingRepair::ResamplingRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
bool ResamplingRepair::resample(Particle* const p, int resamples){
	if (resamples >= 100 || isFeasible(p))
		return false;
	return true;
}

ProjectionBaseRepair::ProjectionBaseRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void ProjectionBaseRepair::repair(Particle* const p, Particle* const base, Particle* const target){
	std::vector<double> x = p->getX();
	std::vector<double>alphas(D+1);
	alphas[D] = 1;

	for (int i = 0; i < D; i++){
		if (x[i] > ub[i]){
			alphas[i] = (base->getX(i) - ub[i])/(base->getX(i) - x[i]);
		} else if (x[i] < lb[i]){
			alphas[i] = (base->getX(i) - lb[i])/(base->getX(i) - x[i]);
		} else
			alphas[i] = std::numeric_limits<double>::max(); //Can't divide by zero
	}

	double alpha=*std::min_element(alphas.begin(), alphas.end());
	if (alpha != 1){
		std::vector<double> b = base->getX();
		scale(b, (1-alpha));
		scale(x, alpha);
		add(x, b, x);
		p->setX(x);
	}
}
