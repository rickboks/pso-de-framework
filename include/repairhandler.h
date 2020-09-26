#pragma once
#include <functional>
#include <vector>
#include "constrainthandler.h"

// Differential Evolution
class RandBaseRepair : public DEConstraintHandler {
	public:
		RandBaseRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){};
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target) const;
};

class MidpointBaseRepair : public DEConstraintHandler {
	public:
		MidpointBaseRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){};
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target) const;
};

class MidpointTargetRepair : public DEConstraintHandler {
	public:
		MidpointTargetRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){};
		void repairDE(Solution* const p, Solution const * const base, Solution const* const target) const;
};

class ConservatismRepair: public DEConstraintHandler {
	public:
		ConservatismRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){};
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target) const;
};

class ProjectionMidpointRepair : public DEConstraintHandler {
	public:
		ProjectionMidpointRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){};
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target) const;
};

class ProjectionBaseRepair: public DEConstraintHandler {
	public:
		ProjectionBaseRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){};
		void repairDE(Solution* const p, Solution const* const base, Solution const* const target) const;
};

class DEReinitializationRepair : public DEConstraintHandler {
	public:
		DEReinitializationRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){}; 
		void repair(Solution* const p) const;
};

class DEProjectionRepair : public DEConstraintHandler {
	public:
		DEProjectionRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){}; 
		void repair(Solution* const p) const;
};

class DEReflectionRepair : public DEConstraintHandler {
	public:
		DEReflectionRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){}; 
		void repair(Solution* const p) const;
};

class DEWrappingRepair : public DEConstraintHandler {
	public:
		DEWrappingRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), DEConstraintHandler(lb, ub){}; 
		void repair(Solution* const p) const;
};

class DETransformationRepair : public DEConstraintHandler { //Adapted from https://github.com/psbiomech/c-cmaes
	private:
		std::vector<double> al, au, xlo, xhi, r;
		void shift(Solution* const p) const;
	public:
		DETransformationRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Solution* const p) const;
};

// Particle Swarm Optimization
class HyperbolicRepair : public PSOConstraintHandler {
	public:
		HyperbolicRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), PSOConstraintHandler(lb, ub){};
		void repairVelocityPre(Particle * const p) const;
};

class PBestDimRepair : public PSOConstraintHandler {
	public:
		PBestDimRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), PSOConstraintHandler(lb, ub){};
		void repairPSO(Particle * const p) const;
};

class PSOReinitializationRepair : public PSOConstraintHandler {
	public:
		PSOReinitializationRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), PSOConstraintHandler(lb, ub){}; 
		void repair(Particle* const p) const;
};

class PSOProjectionRepair : public PSOConstraintHandler {
	public:
		PSOProjectionRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), PSOConstraintHandler(lb, ub){}; 
		void repair(Particle* const p) const;
};

class PSOReflectionRepair : public PSOConstraintHandler {
	public:
		PSOReflectionRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), PSOConstraintHandler(lb, ub){}; 
		void repair(Particle* const p) const;
};

class PSOWrappingRepair : public PSOConstraintHandler {
	public:
		PSOWrappingRepair(std::vector<double> const lb, std::vector<double> const ub):ConstraintHandler(lb,ub), PSOConstraintHandler(lb, ub){}; 
		void repair(Particle* const p) const;
};

class PSOTransformationRepair : public PSOConstraintHandler { //Adapted from https://github.com/psbiomech/c-cmaes
	private:
		std::vector<double> al, au, xlo, xhi, r;
		void shift(Particle* const p) const;
	public:
		PSOTransformationRepair(std::vector<double> const lb, std::vector<double> const ub);
		void repair(Particle* const p) const;
};

// Generic
class ResamplingRepair : public DEConstraintHandler, public PSOConstraintHandler  {
	public:
		ResamplingRepair(std::vector<double> const lb, std::vector<double> const ub)
			:ConstraintHandler(lb,ub), DEConstraintHandler(lb,ub), PSOConstraintHandler(lb, ub){};
		bool resample(Solution const* const p, int const resamples) const;
};

