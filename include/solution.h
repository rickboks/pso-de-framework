#pragma once
#include <vector>
#include <string>
typedef void (*evaluate_function_t)(const double *x, double *y);

class Particle;
class Solution {
	public:
		Solution();
		Solution(int const D);
		Solution(Particle* p);
		Solution(std::vector<double>x); 
		Solution(std::vector<double> x, double fitness);
		//~Solution();
		std::vector<double> getPosition() const;
		double evaluate (evaluate_function_t evalFunc);
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