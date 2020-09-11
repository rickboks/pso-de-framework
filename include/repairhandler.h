#pragma once
#include <vector>
#include "constrainthandler.h"

class Particle;

class RepairHandler : public ConstraintHandler {
	public:
		virtual ~RepairHandler();
		RepairHandler(std::vector<double> const lb, std::vector<double> const ub); 
		virtual void repair(Particle* const p) const; //Generic repair
		virtual void repair(Particle* const p, Particle const* const base, Particle const* const target) const; //DE mutation repair
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

class ResamplingRepair : public RepairHandler {
	public:
		ResamplingRepair(std::vector<double> const lb, std::vector<double> const ub);
		bool resample(Particle const* const p, int const resamples) const;
};
