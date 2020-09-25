#pragma once
#include <vector>
#include <string>
#include <set>
#include <map>
#include <functional>

class Solution;
class Particle;

class ConstraintHandler {
	public:
		ConstraintHandler(std::vector<double> const lb, std::vector<double> const ub);
		virtual ~ConstraintHandler();

		bool isFeasible(Solution const * const p) const;
		virtual bool resample(Solution const * const p, int const resamples) const;

		virtual void repair(Solution* const p) const {};// Generic constraint handler
		virtual void repairDE(Solution* const p, Solution const * const base, Solution const* const target) const {}; // Differential Evolution
		virtual void repairPSO(Particle* const p) const {};// PSO constraint handler
		virtual void repairVelocityPre(Particle* const p) const {};
	protected:
		std::vector<double> const lb;
		std::vector<double> const ub;
		int const D;
};

class DEConstraintHandler : public ConstraintHandler {
	public:
		DEConstraintHandler(std::vector<double>const lb,std::vector<double>const ub);
		virtual ~DEConstraintHandler();
		virtual void repairDE(Solution* const p, Solution const * const base, Solution const* const target) const {}; // Differential Evolution
};

class PSOConstraintHandler : public ConstraintHandler {
	public:
		PSOConstraintHandler(std::vector<double>const lb,std::vector<double>const ub);
		virtual ~PSOConstraintHandler();
		virtual void repairPSO(Particle* const p) const {};// PSO constraint handler
		virtual void repairVelocityPre(Particle* const p) const {};
		//void repairVelocityPost(Solution* const p, int const i) const;
};

extern std::map<std::string, std::function<ConstraintHandler* (std::vector<double>, std::vector<double>)>> const deCHs;
extern std::map<std::string, std::function<ConstraintHandler* (std::vector<double>, std::vector<double>)>> const psoCHs;
