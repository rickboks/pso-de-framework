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
		int nCorrected;
		bool isFeasible(Solution const * const p) const;
	public:
		ConstraintHandler(std::vector<double> const lb, std::vector<double> const ub): lb(lb), ub(ub), D(lb.size()), nCorrected(0){};
		virtual ~ConstraintHandler(){};
		virtual bool resample(Solution* const p, int const resamples);
		virtual void penalize(Solution* const p){};
		int getCorrections() const;
};

class DEConstraintHandler : virtual public ConstraintHandler {
	public:
		DEConstraintHandler(std::vector<double>const lb,std::vector<double>const ub): ConstraintHandler(lb, ub){};
		virtual ~DEConstraintHandler(){};
		virtual void repairDE(Solution* const p, Solution const * const base, Solution const* const target){}; // DE constraint handler
		virtual void repair(Solution* const p){};// Generic constraint handler
};

class PSOConstraintHandler : virtual public ConstraintHandler {
	protected:
		void repairVelocityPost(Particle* const p, int const i); // Change velocity after changing position
	public:
		PSOConstraintHandler(std::vector<double>const lb,std::vector<double>const ub): ConstraintHandler(lb,ub){};
		virtual ~PSOConstraintHandler(){};
		virtual void repairPSO(Particle* const p){};// PSO constraint handler
		virtual void repairVelocityPre(Particle* const p){}; // For constraint handlers that fix the velocity
		virtual void repair(Particle* const p){}; // Generic constraint handler
};

extern std::map<std::string, std::function<DEConstraintHandler* (std::vector<double>, std::vector<double>)>> const deCHs;
extern std::map<std::string, std::function<PSOConstraintHandler* (std::vector<double>, std::vector<double>)>> const psoCHs;
