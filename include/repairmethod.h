#include <vector>
#include "constrainthandler.h"
#include "solution.h"

class RepairHandler: public ContraintHandler {
	protected:
		std::vector<double>const lb;
		std::vector<double>const ub;
		int const D;
	public:
		RepairHandler(std::vector<double>const lb, std::vector<double>const ub);
		virtual void repair(std::vector<double>& x) = 0;
};

class ProjectionRepair : public RepairHandler {
	void repair(std::vector<double>& x) = 0;
	std::vector<Particle*> repair();
};
