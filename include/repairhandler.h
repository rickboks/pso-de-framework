#pragma once
#include <functional>
#include <limits>
#include <vector>
#include "constrainthandler.h"
#include <iostream>

// Differential Evolution
class RandBaseRepair : public DEConstraintHandler {
	public:
		RandBaseRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){};
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target);
};

class MidpointBaseRepair : public DEConstraintHandler {
	public:
		MidpointBaseRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){};
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target);
};

class MidpointTargetRepair : public DEConstraintHandler {
	public:
		MidpointTargetRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){};
		void repairDE(Solution* const p, Solution const * const base, Solution const* const target);
};

class ProjectionMidpointRepair : public DEConstraintHandler {
	public:
		ProjectionMidpointRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){};
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target);
};

class ProjectionBaseRepair: public DEConstraintHandler {
	public:
		ProjectionBaseRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){};
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target);
};

class ConservatismRepair : public DEConstraintHandler {
	public:
		ConservatismRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){};
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target);
};

// Particle Swarm Optimization
class HyperbolicRepair : public PSOConstraintHandler {
	public:
		HyperbolicRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), PSOConstraintHandler(lb, ub){};
		void repairVelocityPre(Particle * const p);
};

class PBestDimRepair : public PSOConstraintHandler {
	public:
		PBestDimRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), PSOConstraintHandler(lb, ub){};
		void repairPSO(Particle * const p);
};

// Generic
class ResamplingRepair : public DEConstraintHandler, public PSOConstraintHandler  {
	public:
		ResamplingRepair(std::vector<double> const lb, std::vector<double> const ub)
			:ConstraintHandler(lb,ub), DEConstraintHandler(lb,ub), PSOConstraintHandler(lb, ub){};
		bool resample(Solution * const p, int const resamples);
};

class DeathPenalty : public DEConstraintHandler, public PSOConstraintHandler {
	public:
		DeathPenalty(std::vector<double>const lb,std::vector<double>const ub):ConstraintHandler(lb,ub), 
			DEConstraintHandler(lb,ub), PSOConstraintHandler(lb,ub){};
		void penalize(Solution* const p);
};

class ReinitializationRepair : public DEConstraintHandler, public PSOConstraintHandler {
	public:
		ReinitializationRepair(std::vector<double> const lb, std::vector<double> const ub)
			:ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub), PSOConstraintHandler(lb,ub){}; 
		void repair(Solution* const p);
		void repair(Particle* const p);
};

class ProjectionRepair : public DEConstraintHandler, public PSOConstraintHandler {
	public:
		ProjectionRepair(std::vector<double> const lb, std::vector<double> const ub)
			:ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub), PSOConstraintHandler(lb,ub){};
		void repair(Solution* const p);
		void repair(Particle* const p);
};

class ReflectionRepair : public DEConstraintHandler, public PSOConstraintHandler {
	public:
		ReflectionRepair(std::vector<double> const lb, std::vector<double> const ub)
			:ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub), PSOConstraintHandler(lb,ub){}; 
		void repair(Solution* const p);
		void repair(Particle* const p);
};

class WrappingRepair : public DEConstraintHandler, public PSOConstraintHandler {
	public:
		WrappingRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub), 
		PSOConstraintHandler(lb,ub){}; 
		void repair(Solution* const p);
		void repair(Particle* const p);
};

class TransformationRepair : public DEConstraintHandler, public PSOConstraintHandler { //Adapted from https://github.com/psbiomech/c-cmaes
	private:
		std::vector<double> al, au, xlo, xhi, r;
		bool shift(Solution* const p);
		bool shift(Particle* const p);
	public:
		TransformationRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Solution* const p);
		void repair(Particle* const p);
};
