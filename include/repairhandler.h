#pragma once
#include <vector>
#include "constrainthandler.h"

class Particle;

class RepairHandler : public PSOConstraintHandler, public DEConstraintHandler {
	protected:
		std::vector<double>const lb;
		std::vector<double>const ub;
		int const D;
	public:
		RepairHandler(std::vector<double>const lb, std::vector<double>const ub);
		virtual void repair(std::vector<double>&x);
		virtual ~RepairHandler()=0;
};

class GenericRepairHandler : public RepairHandler {
	public:
		GenericRepairHandler(std::vector<double> const lb, std::vector<double>const ub);
		virtual ~GenericRepairHandler()=0;
		virtual void repair(std::vector<double>& x)=0;
};

class ReinitializationRepair : public GenericRepairHandler {
	public:
		ReinitializationRepair(std::vector<double> const lb, std::vector<double> const ub); 
	private:
		void repair(std::vector<double>& x);
		std::vector<Particle*> repair();
};

class ProjectionRepair : public GenericRepairHandler {
	public:
		ProjectionRepair(std::vector<double> const lb, std::vector<double> const ub); 
	private:
		void repair(std::vector<double>& x);
		std::vector<Particle*> repair();
};

class ReflectionRepair : public GenericRepairHandler {
	public:
		ReflectionRepair(std::vector<double> const lb, std::vector<double> const ub); 
	private:
		void repair(std::vector<double>& x);
		std::vector<Particle*> repair();
};

class WrappingRepair : public GenericRepairHandler {
	public:
		WrappingRepair(std::vector<double> const lb, std::vector<double> const ub); 
	private:
		void repair(std::vector<double>& x);
		std::vector<Particle*> repair();
};

class ProjectionMidpointRepair : public GenericRepairHandler {
	public:
		ProjectionMidpointRepair(std::vector<double> const lb, std::vector<double> const ub);
	private:
		void repair(std::vector<double>& x);
		std::vector<Particle*> repair();
};
