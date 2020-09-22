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
Particle* getBest(std::vector<Particle*>const& genomes);
Particle* getWorst(std::vector<Particle*>const& genomes);
Particle* pickRandom(std::vector<Particle*>& possibilities);
std::vector<Particle*> pickRandom(std::vector<Particle*>& possibilities, int n);
bool comparePtrs(Particle*a, Particle*b);
void sortOnFitness(std::vector<Particle*>& genomes);
Particle* getPBest(std::vector<Particle*> genomes, double const p);
std::string generateConfig(std::string templateFile, std::string name);
