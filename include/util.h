#pragma once
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include "rng.h"
#include "particle.h"

void scale(std::vector<double>& vec, double x);
void add(std::vector<double>const& lhs, std::vector<double>const& rhs, std::vector<double>& store);
void subtract(std::vector<double>const& lhs, std::vector<double>const& rhs, std::vector<double>& store);
void randomMult(std::vector<double>& vec, double min, double max);
bool comparePtrs(Solution*a, Solution*b);


template <typename T>
void sortOnFitness(std::vector<T*>& genomes){
	std::sort(genomes.begin(), genomes.end(), comparePtrs);
}

template<typename T>
T* getPBest(std::vector<T*> genomes){
	//double const p = rng.randDouble(2./genomes.size(), 0.2);
	double const p = 0.1;
	sortOnFitness(genomes);
	std::vector<T*> bestP = std::vector<T*>(genomes.begin(), genomes.begin() + (genomes.size() * p));
	if (bestP.empty())
		return genomes[0];
	else
		return bestP[rng.randInt(0, bestP.size()-1)];
}

template<typename T>
T* getBest(std::vector<T*>const& genomes){
	int best = 0;
	double bestF = std::numeric_limits<double>::max();

	for (unsigned int i = 0; i < genomes.size(); i++){
		double const score = genomes[i]->getFitness();
		if (score < bestF){
			bestF = score;
			best = i;
		}
	}
	return genomes[best];
}

template<typename T>
T* getWorst(std::vector<T*>const& genomes){
	int worst = 0;
	double worstF = -std::numeric_limits<double>::max();

	for (unsigned int i = 0; i < genomes.size(); i++){
		double const score = genomes[i]->getFitness();
		if (score > worstF){
			worstF = score;
			worst = i;
		}
	}

	return genomes[worst];
}

template<typename T>
T* pickRandom(std::vector<T*>& possibilities){
	int const r = rng.randInt(0,possibilities.size()-1);
	T* g = possibilities[r];
	possibilities.erase(possibilities.begin() + r);
	return g;
}

template<typename T>
std::vector<T*> pickRandom(std::vector<T*>& possibilities, int n){
	std::vector<T*> particles;
	for (int i = 0; i < n; i++){
		particles.push_back(pickRandom(possibilities));
	}
	return particles;
}

std::string generateConfig(std::string templateFile, std::string name);
void printVec(std::vector<double> v);
