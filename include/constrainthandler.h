#pragma once
#include <vector>

class Particle;

class ConstraintHandler {
	protected:
		std::vector<double> const lb;
		std::vector<double> const ub;
		int const D;
	public:
		ConstraintHandler(std::vector<double> lb,std::vector<double> ub);
		virtual ~ConstraintHandler();
		virtual void repair(Particle* const p) const;// Generic constraint handler
		virtual void repair(Particle* const p, Particle const * const base, Particle const* const target) const;
		virtual bool resample(Particle const * const p, int const resamples) const;
		void repairVelocity(Particle* const p, int const i) const;
		bool isFeasible(Particle const * const p) const;
};
