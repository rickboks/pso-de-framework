#pragma once
#include <vector>
#include "constrainthandler.h"

class Particle;

class RepairHandler : public ConstraintHandler {
	public:
		virtual ~RepairHandler();
		RepairHandler(std::vector<double> const lb, std::vector<double> const ub); 
		virtual void repair(Particle* p); //Generic repair
		virtual void repair(Particle* p, Particle* base, Particle* target); //DE mutation repair
		virtual bool resample(Particle* p); //Resampling repair
};   

class ReinitializationRepair : public RepairHandler {
	public:
		ReinitializationRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* p);
};

class ProjectionRepair : public RepairHandler {
	public:
		ProjectionRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* p);
};

class ReflectionRepair : public RepairHandler {
	public:
		ReflectionRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* p);
};

class WrappingRepair : public RepairHandler {
	public:
		WrappingRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* p);
};

// Differential Evolution
class RandBaseRepair : public RepairHandler {
	public:
		RandBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* p, Particle* base, Particle* target);
};

class MidpointBaseRepair : public RepairHandler {
	public:
		MidpointBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* p, Particle* base, Particle* target);
};

class MidpointTargetRepair : public RepairHandler {
	public:
		MidpointTargetRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* p, Particle* base, Particle* target);
};

class ConservatismRepair: public RepairHandler {
	public:
		ConservatismRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* p, Particle* base, Particle* target);
};

class ProjectionMidpointRepair : public RepairHandler {
	public:
		ProjectionMidpointRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* p, Particle* base, Particle* target);
};

class ProjectionBaseRepair: public RepairHandler {
	public:
		ProjectionBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* p, Particle* base, Particle* target);
};

class ResamplingRepair : public RepairHandler {
	public:
		ResamplingRepair(std::vector<double> const lb, std::vector<double> const ub);
		bool resample(Particle* p);
};
