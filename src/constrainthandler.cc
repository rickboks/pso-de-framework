#include "constrainthandler.h"
#include "particle.h"

ConstraintHandler::ConstraintHandler(std::vector<double> lb,std::vector<double> ub)
: lb(lb), ub(ub), D(lb.size()){}

ConstraintHandler::~ConstraintHandler(){}

void ConstraintHandler::repair(Particle* p, Particle* base, Particle* target){}
void ConstraintHandler::repair(Particle* p){}
void ConstraintHandler::repairVelocity(Particle* p, int i){
	if (p->isPSO) p->setV(i, -0.5 * p->getV(i));
}

bool ConstraintHandler::isFeasible(Particle* p){
	for (int i = 0; i < D; i++){
		if (p->getX(i) < lb[i] || p->getX(i) > ub[i]){
			return false;
		}
	}
	return true;
}


//GenericConstraintHandler::GenericConstraintHandler(std::vector<double> lb,std::vector<double> ub)
//: ConstraintHandler(lb,ub){}

//DEConstraintHandler::DEConstraintHandler(std::vector<double> lb,std::vector<double> ub)
//: ConstraintHandler(lb,ub){}

//PSOConstraintHandler::PSOConstraintHandler(std::vector<double> lb,std::vector<double> ub)
//: ConstraintHandler(lb,ub){}


//DEConstraintHandler::~DEConstraintHandler(){}

//PSOConstraintHandler::~PSOConstraintHandler(){}
