#include "rng.h"
#include <algorithm>

RNG::RNG()
: rng(dev()), boolDist(0.5){

} 

bool RNG::randBool(){
	return boolDist(rng);
}

double RNG::randDouble(double start, double end){
	std::uniform_real_distribution<double> dist(start, end);
	return dist(rng);
}

int RNG::randInt(int start, int end){
	std::uniform_int_distribution<int> dist(start, end);
	return dist(rng);
}

double RNG::normalDistribution(double mean, double stdDev){
	std::normal_distribution<double> N(mean, stdDev);
	return N(rng);
}

double RNG::cauchyDistribution(double a, double b){
	std::cauchy_distribution<double> C(a,b);
	return C(rng);
}

void RNG::shuffle(std::vector<double>::iterator first, std::vector<double>::iterator last){
	std::shuffle(first, last, rng);
}

RNG rng; //Global Random Number Generator