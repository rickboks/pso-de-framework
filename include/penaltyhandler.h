#pragma once
#include <vector>
#include "constrainthandler.h"

class Solution;

class DeathPenalty : public DEConstraintHandler, public PSOConstraintHandler {
	public:
		DeathPenalty(std::vector<double>const lb,std::vector<double>const ub);
		void penalize(Solution* const p) const;
};
