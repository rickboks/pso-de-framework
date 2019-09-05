#pragma once
#include <random>
class RNG {
	private:
		std::random_device dev;
		std::mt19937 rng;
		std::bernoulli_distribution boolDist;
	public:
		RNG();
		bool randBool();
		double randDouble(double start, double end);
		int randInt(int start, int end);
		double normalDistribution(double mean, double stdDev);
		double cauchyDistribution(double a, double b);
};

extern RNG rng;