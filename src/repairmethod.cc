#include "repairmethod.h"

RepairHandler::RepairHandler(std::vector<double>const lb, std::vector<double>const ub)
	: lb(lb),ub(ub), D(lb.size()){
}

void ProjectionRepair::repair(std::vector<double>& x){
	for (int i = 0; i < D; i++){
		if (x[i] < lb[i]) 
			x[i] = lb[i];
		else if (x[i] > ub[i]) 
			x[i] = ub[i];
	}
}
