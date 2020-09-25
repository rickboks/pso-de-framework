#pragma once
#include <functional>
#include <vector>
#include "constrainthandler.h"


class Particle;

class RepairHandler : public ConstraintHandler {
	public:
		virtual ~RepairHandler();
		RepairHandler(std::vector<double> const lb, std::vector<double> const ub); 
};   

class ReinitializationRepair : public RepairHandler {
	public:
		ReinitializationRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* const p) const;
};

class ProjectionRepair : public RepairHandler {
	public:
		ProjectionRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* const p) const;
};

class ReflectionRepair : public RepairHandler {
	public:
		ReflectionRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* const p) const;
};

class WrappingRepair : public RepairHandler {
	public:
		WrappingRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* const p) const;
};

class ResamplingRepair : public RepairHandler {
	public:
		ResamplingRepair(std::vector<double> const lb, std::vector<double> const ub);
		bool resample(Particle const* const p, int const resamples) const;
};

class TransformationRepair : public RepairHandler { //Adapted from https://github.com/psbiomech/c-cmaes
	private:
		std::vector<double> al, au, xlo, xhi, r;
		void shift(Particle* const p) const;
	public:
		TransformationRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p) const;
};

// Differential Evolution
class RandBaseRepair : public RepairHandler {
	public:
		RandBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p, Particle const* const base, Particle const* const target) const;
};

class MidpointBaseRepair : public RepairHandler {
	public:
		MidpointBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p, Particle const* const base, Particle const* const target) const;
};

class MidpointTargetRepair : public RepairHandler {
	public:
		MidpointTargetRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p, Particle const * const base, Particle const* const target) const;
};

class ConservatismRepair: public RepairHandler {
	public:
		ConservatismRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p, Particle const* const base, Particle const* const target) const;
};

class ProjectionMidpointRepair : public RepairHandler {
	public:
		ProjectionMidpointRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p, Particle const* const base, Particle const* const target) const;
};

class ProjectionBaseRepair: public RepairHandler {
	public:
		ProjectionBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p, Particle const* const base, Particle const* const target) const;
};

// Particle Swarm Optimization
class HyperbolicRepair : public RepairHandler {
	public:
		HyperbolicRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repairVelocityPre(Particle * const p) const;
};

class PBestDimRepair : public RepairHandler {
	public:
		PBestDimRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle * const p) const;
};
