#pragma once
#include <vector>
class ConstraintHandler {
	public:
	//ConstraintHandler();
	virtual void repair(std::vector<double>& x) = 0;
	virtual ~ConstraintHandler()=0;
};
