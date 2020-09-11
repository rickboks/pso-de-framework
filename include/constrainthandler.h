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
		void repairVelocity(Particle* p, int i);
};

//class DEConstraintHandler : public ConstraintHandler {
	//public:
		//DEConstraintHandler(std::vector<double> lb,std::vector<double> ub);
		//virtual ~DEConstraintHandler();
		//virtual void repair(Particle* p, Particle* base, Particle* target)=0;
		//virtual void repair(Particle* p)=0; //Generic constraint handler
//};

//class PSOConstraintHandler : public ConstraintHandler {
	//public:
		//PSOConstraintHandler(std::vector<double> lb,std::vector<double> ub);
		//virtual ~PSOConstraintHandler();
		//virtual void repair(Particle* p)=0;
//};

//class GenericConstraintHandler : public ConstraintHandler {
	//public:
		//GenericConstraintHandler(std::vector<double> lb,std::vector<double> ub);
		//virtual ~GenericConstraintHandler();
		//virtual void repair(Particle* p)=0;
//};
