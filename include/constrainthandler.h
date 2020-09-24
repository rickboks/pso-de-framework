#pragma once
#include <vector>
#include <string>
#include <set>
#include <map>
#include <functional>

class Particle;
class ConstraintHandler {
	protected:
		std::vector<double> const lb;
		std::vector<double> const ub;
		int const D;
	public:
		ConstraintHandler(std::vector<double>const lb,std::vector<double>const ub);
		virtual ~ConstraintHandler();
		virtual void repair(Particle* const p) const;// Generic constraint handler
		virtual void repair(Particle* const p, Particle const * const base, Particle const* const target) const; // Differential Evolution
		virtual bool resample(Particle const * const p, int const resamples) const;
		virtual void repairVelocityPre(Particle* const p) const;
		void repairVelocityPost(Particle* const p, int const i) const;
		bool isFeasible(Particle const * const p) const;
};

extern std::map<std::string, std::function<ConstraintHandler* (std::vector<double>, std::vector<double>)>> const psoCHs, deCHs;
