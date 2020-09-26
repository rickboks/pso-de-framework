#pragma once
#include <vector>
#include <string>
#include <set>
#include <map>
#include <functional>

class Solution;
class Particle;

class ConstraintHandler {
	protected:
		std::vector<double> const lb;
		std::vector<double> const ub;
		int const D;
		bool isFeasible(Solution const * const p) const;
	public:
		ConstraintHandler(std::vector<double> const lb, std::vector<double> const ub): lb(lb), ub(ub), D(lb.size()){};
		virtual ~ConstraintHandler(){};
		virtual bool resample(Solution const * const p, int const resamples) const;
};

class DEConstraintHandler : virtual public ConstraintHandler {
	public:
		DEConstraintHandler(std::vector<double>const lb,std::vector<double>const ub): ConstraintHandler(lb, ub){};
		virtual ~DEConstraintHandler(){};
		virtual void repairDE(Solution* const p, Solution const * const base, Solution const* const target) const {}; // DE constraint handler
		virtual void repair(Solution* const p) const {};// Generic constraint handler
		virtual void penalize(Solution* const p) const {}; // Generic penalize
};

class PSOConstraintHandler : virtual public ConstraintHandler {
	protected:
		void repairVelocityPost(Particle* const p, int const i) const; // Change velocity after changing position
	public:
		PSOConstraintHandler(std::vector<double>const lb,std::vector<double>const ub): ConstraintHandler(lb,ub){};
		virtual ~PSOConstraintHandler(){};
		virtual void repairPSO(Particle* const p) const {};// PSO constraint handler
		virtual void repairVelocityPre(Particle* const p) const {}; // For constraint handlers that fix the velocity
																	// before updating the position
		virtual void repair(Particle* const p) const {}; // Generic constraint handler
		virtual void penalize(Solution* const p) const {}; // Generic penalize
};

extern std::map<std::string, std::function<DEConstraintHandler* (std::vector<double>, std::vector<double>)>> const deCHs;
extern std::map<std::string, std::function<PSOConstraintHandler* (std::vector<double>, std::vector<double>)>> const psoCHs;
