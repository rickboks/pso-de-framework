#pragma once
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include "rng.h"
#include "particle.h"

void scale(std::vector<double>& vec, double const x);
void add(std::vector<double>const& lhs, std::vector<double>const& rhs, std::vector<double>& store);
void subtract(std::vector<double>const& lhs, std::vector<double>const& rhs, std::vector<double>& store);
void randomMult(std::vector<double>& vec, double const min, double const max);
bool comparePtrs(Solution const* const a, Solution const* const b);
double distance(Solution const*const s1, Solution const*const s2);
std::string generateConfig(std::string const templateFile, std::string const name);
void printVec(std::vector<double> const v);
std::string checkFilename(std::string fn);

template <typename T>
void sortOnFitness(std::vector<T*>& genomes){
	std::sort(genomes.begin(), genomes.end(), comparePtrs);
}

template<typename T>
T* getPBest(std::vector<T*> genomes){
	double const p = std::max(0.05, 3./genomes.size());
	//double const p = 0.1;

	sortOnFitness(genomes);
	std::vector<T*> bestP = std::vector<T*>(genomes.begin(), genomes.begin() + (genomes.size() * p));
	if (bestP.empty())
		return genomes[0];
	else
		return bestP[rng.randInt(0, bestP.size()-1)];
}

template<typename T>
T* getBest(std::vector<T*>const& genomes){
	T* best = NULL;

	for (T* s : genomes)
		if (best == NULL || *s < *best)
			best = s;

	return best;
}

template<typename T>
T* getWorst(std::vector<T*>const& genomes){
	T* worst = NULL;

	for (T* s : genomes)
		if (worst == NULL || *worst < *s)
			worst = s;

	return worst;
}

template<typename T>
T* pickRandom(std::vector<T*>& possibilities){
	int const r = rng.randInt(0,possibilities.size()-1);
	T* g = possibilities[r];
	possibilities.erase(possibilities.begin() + r);
	return g;
}

template<typename T>
std::vector<T*> pickRandom(std::vector<T*>& possibilities, int const n){
	std::vector<T*> particles;
	for (int i = 0; i < n; i++){
		particles.push_back(pickRandom(possibilities));
	}
	return particles;
}

template<typename T>
T* rouletteSelect(std::vector<T*>& possibilities, std::vector<double>& prob){
	double totalProb = std::accumulate(prob.begin(), prob.end(), 0.);
	double rand = rng.randDouble(0.,totalProb);
	for (unsigned int i = 0; i < possibilities.size(); i++){
		rand -= prob[i];
		if (rand <= 0.){
			T* selected = possibilities[i];
			possibilities.erase(possibilities.begin() + i);
			prob.erase(prob.begin() + i);
			return selected;
		} 
	}
	return NULL; //should not happen
}

template<typename T>
std::vector<T*> rouletteSelect(std::vector<T*>& possibilities, std::vector<double>& prob, int const n){
	std::vector<T*> particles;
	for (int i = 0; i < n; i++){
		particles.push_back(rouletteSelect(possibilities, prob));
	}
	return particles;
}
