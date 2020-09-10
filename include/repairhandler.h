#pragma once
#include <vector>
#include "constrainthandler.h"

class Particle;

class RepairHandler {
	protected:
		std::vector<double>const lb;
		std::vector<double>const ub;
		int const D;
	public:
		RepairHandler(std::vector<double>const lb, std::vector<double>const ub);
		virtual void repair(Particle* p);
		virtual void repair(Particle* p, Particle* base, Particle* target);
		virtual ~RepairHandler()=0;
		void repairVelocity(Particle* p, int i);
};

// Generic
class GenericRepairHandler : public RepairHandler, public DEConstraintHandler, public PSOConstraintHandler {
	public:
		GenericRepairHandler(std::vector<double> const lb, std::vector<double>const ub);
		virtual ~GenericRepairHandler()=0;
		virtual void repair(Particle* p)=0; // Generic
		void repair(Particle* p, Particle* base, Particle* target); //DE
};

class ReinitializationRepair : public GenericRepairHandler {
	public:
		ReinitializationRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* p);
};

class ProjectionRepair : public GenericRepairHandler {
	public:
		ProjectionRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* p);
};

class ReflectionRepair : public GenericRepairHandler {
	public:
		ReflectionRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* p);
};

class WrappingRepair : public GenericRepairHandler {
	public:
		WrappingRepair(std::vector<double> const lb, std::vector<double> const ub); 
		void repair(Particle* p);
};

// Differential Evolution
class DERepairHandler : public RepairHandler, public DEConstraintHandler {
	public:
		DERepairHandler(std::vector<double> const lb, std::vector<double>const ub);
		virtual ~DERepairHandler()=0;
		virtual void repair(Particle* p, Particle* base, Particle* target)=0;
		void repair(Particle* p);
};

class RandBaseRepair : public DERepairHandler {
	public:
		RandBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* p, Particle* base, Particle* target);
};

class MidpointBaseRepair : public DERepairHandler {
	public:
		MidpointBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* p, Particle* base, Particle* target);
};

class MidpointTargetRepair : public DERepairHandler {
	public:
		MidpointTargetRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* p, Particle* base, Particle* target);
};

class ConservatismRepair: public DERepairHandler {
	public:
		ConservatismRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* p, Particle* base, Particle* target);
};

class ProjectionMidpointRepair : public DERepairHandler {
	public:
		ProjectionMidpointRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* p);
};

class ProjectionBaseRepair: public DERepairHandler {
	public:
		ProjectionBaseRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* p, Particle* base, Particle* target);
};
