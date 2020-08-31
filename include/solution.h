#pragma once
#include <vector>
#include <string>
#include <memory>

class Particle;
template <typename T>
class IOHprofiler_problem;
class IOHprofiler_csv_logger;

class Solution {
	public:
		Solution();
		Solution(int const D);
		Solution(Particle* p);
		Solution(std::vector<double>x); 
		Solution(std::vector<double> x, double fitness);
		std::vector<double> getPosition() const;
		double evaluate (std::shared_ptr<IOHprofiler_problem<double> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger);
		double getFitness() const;
		std::string positionString() const;
		void setPosition(std::vector<double> position, double fitness);
		void randomize(std::vector<double> lowerBounds, std::vector<double> upperBounds);
		int getDimension();
		bool operator < (const Solution& s) const;
	protected:
		std::vector<double> x;
		int const D;
		bool evaluated;
		double fitness;
};
