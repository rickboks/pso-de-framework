#include "genome.h"
#include "rng.h"
#include <random>
#include <algorithm>
#include <functional>
#include "particle.h"
#include <iostream>

Genome::Genome(int const D): 
	Solution(D){
}

Genome::Genome(Particle* particle) :
	Solution(particle){
}

Genome::Genome(std::vector<double>x): 
	Solution(x){
}
