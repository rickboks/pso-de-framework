#pragma once
#include <vector>

class Particle;

class ConstraintHandler {
	protected:
		std::vector<double>const lb;
		std::vector<double>const ub;
		int const D;
	public:
		ConstraintHandler(std::vector<double> lb,std::vector<double> ub);
		virtual ~ConstraintHandler();
		virtual void repair(Particle* p)=0;// Generic constraint handler
		virtual void repair(Particle* p, Particle* base, Particle* target)=0;
		virtual bool resample(Particle* p, int resamples)=0;
		void repairVelocity(Particle* p, int i);
		bool isFeasible(Particle* p);
};
