#include "repairhandler.h"
#include "particle.h"
#include "util.h"
#include <algorithm>
#include <limits>

RepairHandler::~RepairHandler(){}
RepairHandler::RepairHandler(std::vector<double>const lb, std::vector<double>const ub) : ConstraintHandler(lb, ub){}

ReinitializationRepair::ReinitializationRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void ReinitializationRepair::repair(Particle* const p) const{
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i] || p->getX(i) > ub[i]){
			p->setX(i, rng.randDouble(lb[i], ub[i]));
			repairVelocityPost(p, i);
		}
	}
}

ProjectionRepair::ProjectionRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void ProjectionRepair::repair(Particle* const p) const{
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i]){
			p->setX(i, lb[i]);
			repairVelocityPost(p, i);
		} else if (p->getX(i) > ub[i]){
			p->setX(i, ub[i]);
			repairVelocityPost(p, i);
		}
	}
}

ReflectionRepair::ReflectionRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void ReflectionRepair::repair(Particle* const p) const{
	for (int i = 0; i < D; i++){
		bool is_repaired = false;
		while (p->getX(i) < lb[i]){
			p->setX(i, 2. * lb[i] - p->getX(i));
			is_repaired = true;
		}
		while (p->getX(i) > ub[i]){
			p->setX(i, 2. * ub[i] - p->getX(i));
			is_repaired = true;
		}
		if (is_repaired)
			repairVelocityPost(p, i);
	}
}

WrappingRepair::WrappingRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void WrappingRepair::repair(Particle* const p) const{
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i]){
			p->setX(i, ub[i] - std::fmod(lb[i] - p->getX(i), std::abs(ub[i]-lb[i])));
			repairVelocityPost(p,i);
		} else if (p->getX(i) > ub[i]) {
			p->setX(i, lb[i] + std::fmod(p->getX(i) - ub[i], std::abs(ub[i]-lb[i])));
			repairVelocityPost(p,i);
		}
	}
}

ResamplingRepair::ResamplingRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
bool ResamplingRepair::resample(Particle const* const p, int const resamples) const{
	if (resamples >= 100 || isFeasible(p))
		return false;
	return true;
}

//Adapted from https://github.com/psbiomech/c-cmaes
TransformationRepair::TransformationRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub), 
	al(D), au(D), xlo(D), xhi(D), r(D){
	for (int i = 0; i < D; i++){
		al[i] = std::min( (ub[i]-lb[i])/2., (1.+std::abs(lb[i]))/20. );
		au[i] = std::min( (ub[i]-lb[i])/2., (1.+std::abs(ub[i]))/20. );
		xlo[i] = lb[i] - 2. * al[i] - (ub[i] - lb[i]) / 2.;
		xhi[i] = ub[i] + 2. * au[i] + (ub[i] - lb[i]) / 2.;
		r[i] = 2.*(ub[i] - lb[i] + al[i] + au[i]);
	}
}
void TransformationRepair::repair(Particle* const p) const{
	shift(p);
	for (int i = 0; i < D; i++){
		double const x_i = p->getX(i);
		if (x_i < lb[i] + al[i])
			p->setX(i, lb[i] + pow(x_i - (lb[i] - al[i]),2.)/(4.*al[i]));
		else if (x_i > ub[i]-au[i])
			p->setX(i, ub[i] - pow(x_i - (ub[i] + au[i]),2.)/(4.*au[i]));
	}
}
void TransformationRepair::shift(Particle* const p) const {
	for (int i = 0; i < D; i++){
		if (p->getX(i) < xlo[i]) // Shift up
			p->setX(i, p->getX(i) + r[i] * (1 + (int)((xlo[i] - p->getX(i))/r[i])) );
		if (p->getX(i) > xhi[i]) // Shift down
			p->setX(i, p->getX(i) - r[i] * (1 + (int)((p->getX(i)-xhi[i])/r[i])) );

		if (p->getX(i) < lb[i] - al[i]) // mirror
			p->setX(i, p->getX(i) + 2. * (lb[i] - al[i] - p->getX(i)));
		if (p->getX(i) > ub[i] + au[i])
			p->setX(i, p->getX(i) - 2. * (p->getX(i) - ub[i] - au[i]));
	}
}

// Differential Evolution
RandBaseRepair::RandBaseRepair(std::vector<double>const lb, std::vector<double>const ub): RepairHandler(lb, ub){}
void RandBaseRepair::repair(Particle* const p, Particle const* const base, Particle const* const target) const{
	for (int i = 0; i < D; i++){
		if (p->getX(i) > ub[i]){
			p->setX(i, base->getX(i) + rng.randDouble(0,1) * (ub[i] - base->getX(i)));
		} else if (p->getX(i) < lb[i]){
			p->setX(i, base->getX(i) + rng.randDouble(0,1) * (lb[i] - base->getX(i)));
		}
	}
}

MidpointBaseRepair::MidpointBaseRepair(std::vector<double>const lb, std::vector<double>const ub): RepairHandler(lb, ub){}
void MidpointBaseRepair::repair(Particle* const p, Particle const* const base, Particle const* const target) const{
	for (int i = 0; i < D; i++){
		if (p->getX(i) > ub[i]){
			p->setX(i, 0.5 * (base->getX(i) + ub[i]));
		} else if (p->getX(i) < lb[i]){
			p->setX(i, 0.5 * (base->getX(i) + lb[i]));
		}
	}
}

MidpointTargetRepair::MidpointTargetRepair(std::vector<double>const lb, std::vector<double>const ub): RepairHandler(lb, ub){}
void MidpointTargetRepair::repair(Particle* const p, Particle const* const base, Particle const* const target) const{
	for (int i = 0; i < D; i++){
		if (p->getX(i) > ub[i]){
			p->setX(i, 0.5 * (target->getX(i) + ub[i]));
		} else if (p->getX(i) < lb[i]){
			p->setX(i, 0.5 * (target->getX(i) + lb[i]));
		}
	}
}

ConservatismRepair::ConservatismRepair(std::vector<double>const lb, std::vector<double>const ub): RepairHandler(lb, ub){}
void ConservatismRepair::repair(Particle* const p, Particle const* const base, Particle const* const target) const{
	for (int i = 0; i < D; i++){
		if (p->getX(i) > ub[i] || p->getX(i) < lb[i]){
			p->setX(i, base->getX(i));
		}
	}
}

ProjectionMidpointRepair::ProjectionMidpointRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void ProjectionMidpointRepair::repair(Particle* const p, Particle const* const base, Particle const* const target) const{
	std::vector<double> x = p->getX();
	std::vector<double>alphas(D+1);
	alphas[D] = 1;

	for (int i = 0; i < D; i++){
		if (x[i] > ub[i]){
			alphas[i] = (lb[i] - ub[i])/(lb[i] - 2. * x[i] + ub[i]);
		} else if (x[i] < lb[i]){
			alphas[i] = (ub[i] - lb[i])/(lb[i] - 2. * x[i] + ub[i]);
		} else
			alphas[i] = std::numeric_limits<double>::max(); //Can't divide by zero
	}

	double alpha=*std::min_element(alphas.begin(), alphas.end());
	if (alpha != 1.){
		std::vector<double> midpoint(D);
		add(lb, ub, midpoint);
		scale(midpoint, 0.5*(1.-alpha));
		scale(x, alpha);
		add(x, midpoint, x);
		p->setX(x);
	}
}

ProjectionBaseRepair::ProjectionBaseRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void ProjectionBaseRepair::repair(Particle* const p, Particle const* const base, Particle const* const target) const{
	std::vector<double> x = p->getX();
	std::vector<double> alphas(D+1);
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
	if (alpha != 1.){
		std::vector<double> b = base->getX();
		scale(b, (1.-alpha));
		scale(x, alpha);
		add(x, b, x);
		p->setX(x);
	}
}

// Particle Swarm Optimization
HyperbolicRepair::HyperbolicRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void HyperbolicRepair::repairVelocityPre(Particle * const p) const{
	for (int i = 0; i < D; i++){
		double const center = (lb[i] + ub[i])/2.;
		double const v = p->getV(i);
		if (v > center){
			p->setV(i, v / (1. + std::abs(v / (ub[i] - p->getX(i)))));
		} else {
			p->setV(i, v / (1. + std::abs(v / (p->getX(i) - lb[i]))));
		}
	}
}

PBestDimRepair::PBestDimRepair(std::vector<double>const lb, std::vector<double>const ub) : RepairHandler(lb, ub){}
void PBestDimRepair::repair(Particle* const p) const {
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i] || p->getX(i) > ub[i]){
			p->setX(i, p->getP(i));
			p->setV(i,0.);
		}
	}	
}
