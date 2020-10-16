#pragma once
#include <vector>
#include <IOHprofiler_experimenter.h>

class Solution {
	protected:
		std::vector<double> x;
		bool evaluated;
		double fitness;
	public:
		Solution(int const D);
		virtual ~Solution();
		Solution(std::vector<double> const mutant);
		int const D;
		virtual void setX(std::vector<double> const x, double const fitness);
		void setX(int const dim, double const val);
		void setX(std::vector<double> x);
		std::vector<double> getX() const;
		double getX(int const dim) const;
		double evaluate (std::shared_ptr<IOHprofiler_problem<double> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger);
		double getFitness() const;
		void setFitness(double const d);
		std::string positionString() const;
		void randomize(std::vector<double> const lowerBounds, std::vector<double> const upperBounds);
		bool operator < (Solution const& s) const;
		void copy (Solution const * const other);
};
