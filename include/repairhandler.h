#pragma once
#include <vector>
#include "constrainthandler.h"

class Particle;

class RepairHandler : public ConstraintHandler {
	public:
		virtual ~RepairHandler();
		RepairHandler(std::vector<double> const lb, std::vector<double> const ub); 
		virtual void repair(Particle* const p); //Generic repair
		virtual void repair(Particle* const p, Particle* const base, Particle* const target); //DE mutation repair
};   

class ReinitializationRepair : public RepairHandler {
	public:
		ReinitializationRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* const p);
};

class ProjectionRepair : public RepairHandler {
	public:
		ProjectionRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* const p);
};

class ReflectionRepair : public RepairHandler {
	public:
		ReflectionRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* const p);
};

class WrappingRepair : public RepairHandler {
	public:
		WrappingRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* const p);
};

// Differential Evolution
class RandBaseRepair : public RepairHandler {
	public:
		RandBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p, Particle* const base, Particle* const target);
};

class MidpointBaseRepair : public RepairHandler {
	public:
		MidpointBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p, Particle* const base, Particle* const target);
};

class MidpointTargetRepair : public RepairHandler {
	public:
		MidpointTargetRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p, Particle* const base, Particle* const target);
};

class ConservatismRepair: public RepairHandler {
	public:
		ConservatismRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p, Particle* const base, Particle* const target);
};

class ProjectionMidpointRepair : public RepairHandler {
	public:
		ProjectionMidpointRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p, Particle* const base, Particle* const target);
};

class ProjectionBaseRepair: public RepairHandler {
	public:
		ProjectionBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p, Particle* const base, Particle* const target);
};

class ResamplingRepair : public RepairHandler {
	public:
		ResamplingRepair(std::vector<double> const lb, std::vector<double> const ub);
		bool resample(Particle* const p, int const resamples);
};
