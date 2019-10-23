#pragma once
#include <random>
#include <algorithm>
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
		void shuffle(std::vector<double>::iterator first, std::vector<double>::iterator last);
};

extern RNG rng;