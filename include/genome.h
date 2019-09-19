#pragma once
#include <vector>
#include "coco.h"
typedef void (*evaluate_function_t)(const double *x, double *y);
#include "crossovermanager.h"

class Particle;
class Genome {
	private:
		int const D;
		std::vector<double> x;
		bool evaluated;
		double fitness;
		double weightFactor; //used only for population-based mutation
	public:
		Genome(int const D);
		Genome(Particle* particle);
		Genome(std::vector<double>x);
		Genome(double const Cr, CrossoverType const crossoverType, Genome const* parent, Genome const* mutant);
		void randomize(std::vector<double> lowerBounds, std::vector<double> upperBounds);
		std::vector<double> getPosition() const;
		void setPosition(std::vector<double> const x);
		void print() const;
		int getDimension() const;
		double getFitness() const;
		double evaluate(evaluate_function_t evalFunc);

		void setWeightFactor(double weightFactor);
		double getWeightFactor() const;

		bool operator < (const Genome& genome) const;
};