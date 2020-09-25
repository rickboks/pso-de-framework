#pragma once
#include <functional>
#include <vector>
#include "constrainthandler.h"


class Particle;

class ReinitializationRepair : public ConstraintHandler {
	public:
		ReinitializationRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Solution* const p) const;
};

class ProjectionRepair : public ConstraintHandler {
	public:
		ProjectionRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Solution* const p) const;
};

class ReflectionRepair : public ConstraintHandler {
	public:
		ReflectionRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Solution* const p) const;
};

class WrappingRepair : public ConstraintHandler {
	public:
		WrappingRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Solution* const p) const;
};

class ResamplingRepair : public ConstraintHandler {
	public:
		ResamplingRepair(std::vector<double> const lb, std::vector<double> const ub);
		bool resample(Solution const* const p, int const resamples) const;
};

class TransformationRepair : public ConstraintHandler { //Adapted from https://github.com/psbiomech/c-cmaes
	private:
		std::vector<double> al, au, xlo, xhi, r;
		void shift(Solution* const p) const;
	public:
		TransformationRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Solution* const p) const;
};

// Differential Evolution
class RandBaseRepair : public DEConstraintHandler {
	public:
		RandBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target) const;
};

class MidpointBaseRepair : public DEConstraintHandler {
	public:
		MidpointBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target) const;
};

class MidpointTargetRepair : public DEConstraintHandler {
	public:
		MidpointTargetRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repairDE(Solution* const p, Solution const * const base, Solution const* const target) const;
};

class ConservatismRepair: public DEConstraintHandler {
	public:
		ConservatismRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target) const;
};

class ProjectionMidpointRepair : public DEConstraintHandler {
	public:
		ProjectionMidpointRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target) const;
};

class ProjectionBaseRepair: public DEConstraintHandler {
	public:
		ProjectionBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target) const;
};

// Particle Swarm Optimization
class HyperbolicRepair : public PSOConstraintHandler {
	public:
		HyperbolicRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repairVelocityPre(Particle * const p) const;
};

class PBestDimRepair : public PSOConstraintHandler {
	public:
		PBestDimRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repairPSO(Particle * const p) const;
};
