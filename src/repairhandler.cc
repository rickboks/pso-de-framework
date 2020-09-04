#include "repairhandler.h"
#include "util.h"
#include <algorithm>
#include <limits>

RepairHandler::RepairHandler(std::vector<double>const lb, std::vector<double>const ub)
	: PSOConstraintHandler(), DEConstraintHandler(), lb(lb),ub(ub), D(lb.size()){
}

RepairHandler::~RepairHandler(){}

void RepairHandler::repair(std::vector<double>& x){
	// do nothing
}

GenericRepairHandler::GenericRepairHandler(std::vector<double>const lb, std::vector<double>const ub): RepairHandler(lb,ub){}

GenericRepairHandler::~GenericRepairHandler(){}

ReinitializationRepair::ReinitializationRepair(std::vector<double>const lb, std::vector<double>const ub)
	: GenericRepairHandler(lb, ub){
}
void ReinitializationRepair::repair(std::vector<double>& x){
	for (int i = 0; i < D; i++){
		if (x[i] < lb[i] || x[i] > ub[i])
			x[i] = rng.randDouble(lb[i], ub[i]);
	}
}

ProjectionRepair::ProjectionRepair(std::vector<double>const lb, std::vector<double>const ub)
	: GenericRepairHandler(lb, ub){
}

void ProjectionRepair::repair(std::vector<double>& x){
	for (int i = 0; i < D; i++){
		if (x[i] < lb[i]) 
			x[i] = lb[i];
		else if (x[i] > ub[i]) 
			x[i] = ub[i];
	}
}

ReflectionRepair::ReflectionRepair(std::vector<double>const lb, std::vector<double>const ub)
	: GenericRepairHandler(lb, ub){
}
void ReflectionRepair::repair(std::vector<double>& x){
	for (int i = 0; i < D; i++){
		while (true){
			if (x[i] < lb[i]) 
				x[i] = 2 * lb[i] - x[i];
			else if (x[i] > ub[i]) 
				x[i] = 2 * ub[i] - x[i];
			else 
				break;
		}
	}
}

WrappingRepair::WrappingRepair(std::vector<double>const lb, std::vector<double>const ub)
	: GenericRepairHandler(lb, ub){
}
void WrappingRepair::repair(std::vector<double>& x){
	for (int i = 0; i < D; i++){
		while (true){
			if (x[i] < lb[i]) 
				x[i] = ub[i] + x[i] - lb[i];
			else if (x[i] > ub[i]) 
				x[i] = lb[i] + x[i] - ub[i];
			else 
				break;
		}
	}
}

/* Projection to Midpoint */
ProjectionMidpointRepair::ProjectionMidpointRepair(std::vector<double>const lb, std::vector<double>const ub)
	: GenericRepairHandler(lb, ub){
}

void ProjectionMidpointRepair::repair(std::vector<double>& x){
	std::vector<double>alphas(D+1);
	alphas[D] = 1;

	for (int i = 0; i < D; i++){
		if (x[i] > 0)
			alphas[i] = ub[i]/x[i];
		else if (x[i] < 0)
			alphas[i] = lb[i]/x[i];
		else
			alphas[i] = std::numeric_limits<double>::max(); //Can't divide by zero
	}

	double alpha=*std::min_element(alphas.begin(), alphas.end());
	scale(x, alpha);
}
