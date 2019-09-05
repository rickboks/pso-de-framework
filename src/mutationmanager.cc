#include "mutationmanager.h"
#include <algorithm>
#include <iostream>
#include <functional>
#include <limits>
#include "rng.h"
#include "vectoroperations.h"

// FACTORY
MutationManager* MutationManagerFactory::createMutationManager(MutationType const mutationType, std::vector<Genome*>& genomes, double const F){
	switch(mutationType){
		case RAND_1:
			return new Rand1MutationManager(genomes,F);
		case TTB_1:
			return new TTB1MutationManager(genomes, F);
		case BEST_1:
			return new Best1MutationManager(genomes, F);
		case BEST_2:
			return new Best2MutationManager(genomes,F);
		case RAND_2:
			return new Rand2MutationManager(genomes,F);
		case RAND_2_DIR:
			return new Rand2DirMutationManager(genomes,F);
		case NSDE:
			return new NSDEMutationManager(genomes,F);
		case TOPOLOGY:
			return new TopologyMutationManager(genomes,F);

		default:
			throw std::invalid_argument("Error: Invalid DE mutation type");
	}
}

//BASE
MutationManager::MutationManager(std::vector<Genome*>& genomes, double const F)
:genomes(genomes), F(F), D(genomes[0]->getDimension()), popSize(genomes.size()){

}

MutationManager::~MutationManager(){

}

Genome* MutationManager::getBest(){
	int best = 0;
	double bestF = std::numeric_limits<double>::max();

	for (int i = 0; i < popSize; i++){
		double const score = genomes[i]->getFitness();
		if (score < bestF){
			bestF = score;
			best = i;
		}
	}

	return genomes[best];
}

Genome* MutationManager::pickRandom(std::vector<Genome*>& possibilities) {
	int const r = rng.randInt(0,possibilities.size()-1);
	Genome* g = possibilities[r];
	possibilities.erase(possibilities.begin() + r);
	return g;
}

//RAND1
Rand1MutationManager::Rand1MutationManager(std::vector<Genome*>& genomes, double const F)
:MutationManager(genomes, F){

}

std::vector<Genome*> Rand1MutationManager::mutate(){
	std::vector<Genome*> mutants;
	mutants.reserve(popSize);

	for (int i = 0; i < popSize; i++){
		std::vector<Genome*> possibilities = genomes;

		auto position = std::find(possibilities.begin(), possibilities.end(), genomes[i]);
		if (position != possibilities.end())
		    possibilities.erase(position);

		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();
		std::vector<double> xr2 = pickRandom(possibilities)->getX();

		Genome* mutant = new Genome(D);

		std::vector<double> x = xr0;
		std::vector<double> subtraction(D);
		subtract(xr1, xr2, subtraction);
		multiply(subtraction, F);
		add(x,subtraction, x);

		mutant->setX(x);
		mutants.push_back(mutant);		
	}
	return mutants;
}

// TTB1
TTB1MutationManager::TTB1MutationManager(std::vector<Genome*>& genomes, double const F)
:MutationManager(genomes, F){

}

std::vector<Genome*> TTB1MutationManager::mutate(){
	std::vector<double> best = getBest()->getX(); 
	std::vector<Genome*> mutants;
	mutants.reserve(popSize);

	for (int i = 0; i < popSize; i++){
		std::vector<Genome*> possibilities = genomes;

		auto position = std::find(possibilities.begin(), possibilities.end(), genomes[i]);
		if (position != possibilities.end())
		    possibilities.erase(position);

		std::vector<double> mutant = genomes[i]->getX();
		std::vector<double> subtraction(D);
		subtract(best, genomes[i]->getX(), subtraction);
		multiply(subtraction, F);
		add(mutant, subtraction, mutant);

		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();

		subtract(xr0, xr1, subtraction);
		multiply(subtraction,F);
		add(mutant, subtraction, mutant);

		Genome* m = new Genome(D);
		m->setX(mutant);

		mutants.push_back(m);
	}
	
	return mutants;
}


// BEST 1
Best1MutationManager::Best1MutationManager(std::vector<Genome*>& genomes, double const F)
	:MutationManager(genomes, F){

}

std::vector<Genome*> Best1MutationManager::mutate(){
	std::vector<Genome*> mutants;
	mutants.reserve(popSize);
	std::vector<double> best = getBest()->getX();

	for (int i = 0; i < popSize; i++){
		std::vector<Genome*> possibilities = genomes;

		auto position = std::find(possibilities.begin(), possibilities.end(), genomes[i]);
		if (position != possibilities.end())
		    possibilities.erase(position);

		std::vector<double> mutant = best;
		std::vector<double> subtraction(D);

		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();

		subtract(xr0, xr1, subtraction);
		multiply(subtraction,F);
		add(mutant, subtraction, mutant);

		Genome* m = new Genome(D);
		m->setX(mutant);

		mutants.push_back(m);
	}

	return mutants;
}

// BEST 2
Best2MutationManager::Best2MutationManager(std::vector<Genome*>& genomes, double const F)
	:MutationManager(genomes, F){

}

std::vector<Genome*> Best2MutationManager::mutate(){
	std::vector<double> best = getBest()->getX(); 
	std::vector<Genome*> mutants;
	mutants.reserve(popSize);

	for (int i = 0; i < popSize; i++){
		std::vector<Genome*> possibilities = genomes;

		auto position = std::find(possibilities.begin(), possibilities.end(), genomes[i]);
		if (position != possibilities.end())
		    possibilities.erase(position);

		std::vector<double> mutant = best;
		std::vector<double> subtraction(D);

		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();
		std::vector<double> xr2 = pickRandom(possibilities)->getX();
		std::vector<double> xr3 = pickRandom(possibilities)->getX();
		subtract(xr0, xr1, subtraction);
		multiply(subtraction, F);
		add(mutant, subtraction, mutant);
		subtract(xr2, xr3, subtraction);
		multiply(subtraction, F);
		add(mutant, subtraction, mutant);

		Genome* m = new Genome(D);
		m->setX(mutant);

		mutants.push_back(m);
	}

	return mutants;
}

// RAND 2
Rand2MutationManager::Rand2MutationManager(std::vector<Genome*>& genomes, double const F)
	:MutationManager(genomes, F){

}

std::vector<Genome*> Rand2MutationManager::mutate(){
	std::vector<Genome*> mutants;
	mutants.reserve(popSize);

	for (int i = 0; i < popSize; i++){
		std::vector<Genome*> possibilities = genomes;

		auto position = std::find(possibilities.begin(), possibilities.end(), genomes[i]);
		if (position != possibilities.end())
		    possibilities.erase(position);		

		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();
		std::vector<double> xr2 = pickRandom(possibilities)->getX();
		std::vector<double> xr3 = pickRandom(possibilities)->getX();
		std::vector<double> xr4 = pickRandom(possibilities)->getX();

		std::vector<double> mutant = xr4;
		std::vector<double> subtraction(D);
		subtract(xr0, xr1, subtraction);
		multiply(subtraction, F);
		add(mutant, subtraction, mutant);
		subtract(xr2, xr3, subtraction);
		multiply(subtraction, F);
		add(mutant, subtraction, mutant);

		Genome* m = new Genome(D);
		m->setX(mutant);

		mutants.push_back(m);
	}
	
	return mutants;
}

// RAND 2 DIR
Rand2DirMutationManager::Rand2DirMutationManager(std::vector<Genome*>& genomes, double const F)
	:MutationManager(genomes, F){

}

std::vector<Genome*> Rand2DirMutationManager::mutate(){
	std::vector<Genome*> mutants;
	mutants.reserve(popSize);

	for (int i = 0; i < popSize; i++){
		std::vector<Genome*> possibilities = genomes;

		auto position = std::find(possibilities.begin(), possibilities.end(), genomes[i]);
		if (position != possibilities.end())
		    possibilities.erase(position);		

		Genome* r0 = pickRandom(possibilities);
		Genome* r1 = pickRandom(possibilities);
		Genome* r2 = pickRandom(possibilities);
		Genome* r3 = pickRandom(possibilities);

		std::vector<double> xr0, xr1, xr2, xr3;

		if (r0->getFitness() < r1->getFitness()){
			xr0 = r0->getX();
			xr1 = r1->getX();
		} else {
			xr0 = r1->getX();
			xr1 = r0->getX();
		}

		if (r2->getFitness() < r3->getFitness()){
			xr2 = r2->getX();
			xr3 = r3->getX();
		} else {
			xr2 = r3->getX();
			xr3 = r2->getX();
		}

		std::vector<double> mutant = xr0;
		std::vector<double> subtraction(D);
		subtract(xr0, xr1, subtraction);
		add(subtraction, xr2, subtraction);
		subtract(subtraction, xr3, subtraction);
		multiply(subtraction, F*0.5);
		add(mutant, subtraction, mutant);

		Genome* m = new Genome(D);
		m->setX(mutant);

		mutants.push_back(m);
	}
	
	return mutants;
}

// NSDE

NSDEMutationManager::NSDEMutationManager(std::vector<Genome*>& genomes, double const F)
	:MutationManager(genomes, F){

}

std::vector<Genome*> NSDEMutationManager::mutate(){
	std::vector<Genome*> mutants;
	mutants.reserve(popSize);

	for (int i = 0; i < popSize; i++){
		std::vector<Genome*> possibilities = genomes;

		auto position = std::find(possibilities.begin(), possibilities.end(), genomes[i]);
		if (position != possibilities.end())
		    possibilities.erase(position);

		std::vector<double> mutant = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();
		std::vector<double> xr2 = pickRandom(possibilities)->getX();

		std::vector<double> subtraction(D);
		subtract(xr1, xr2, subtraction);

		if (rng.randDouble(0,1) < 0.5){
			double n = rng.normalDistribution(0.5,0.5);
			multiply(subtraction, n);
		} else {
			// std::cauchy_distribution<double> cauchy(0,1);
			double c = rng.cauchyDistribution(0,1);
			multiply(subtraction, c);
		}

		add(mutant, subtraction, mutant);

		Genome* m = new Genome(D);
		m->setX(mutant);

		mutants.push_back(m);
	}

	return mutants;
}

// TOPOLOGY BASED

std::vector<Genome*> TopologyMutationManager::getNeighbors(int const i) const{
	std::vector<Genome*> neighbors;
	neighbors.reserve(2*radius);

	for (int k = 1; k <= radius; k++){
		int index = i - k;
		if (index < 0)
			index = index + popSize;
		neighbors.push_back(genomes[index]);
	}

	for (int k = 1; k <= radius; k++){
		int index = (i + k) % popSize;
		neighbors.push_back(genomes[index]);
	}
	return neighbors;
}

TopologyMutationManager::TopologyMutationManager(std::vector<Genome*>& genomes, double const F)
	:MutationManager(genomes, F), radius(3), alpha(F), beta(F){
	for (Genome* g : genomes){
		g->setWeightFactor(rng.randDouble(0,1));
	}
}

std::vector<Genome*> TopologyMutationManager::mutate(){
	std::vector<Genome*> mutants;
	mutants.reserve(popSize);

	Genome* best = getBest();
	std::vector<double> bestX = best->getX();

	for (int i = 0; i < (int)genomes.size(); i++){
		std::vector<Genome*> neighbors = getNeighbors(i);

		double bestF = std::numeric_limits<double>::max();

		std::vector<double> bestNeighbor;
		for (Genome* g : neighbors){
			if (g->getFitness() < bestF){
				bestF = g->getFitness();
				bestNeighbor = g->getX();
			}
		}

		// LOCAL VECTOR CREATION
		std::vector<double> xr0 = pickRandom(neighbors)->getX();
		std::vector<double> xr1 = pickRandom(neighbors)->getX();

		std::vector<double> localVector = genomes[i]->getX();
		std::vector<double> subtraction(D);
		subtract(bestNeighbor, genomes[i]->getX(), subtraction);
		multiply(subtraction, alpha);
		add(subtraction, localVector, localVector);
		subtract(xr0, xr1, subtraction);
		multiply(subtraction, beta);
		add(subtraction, localVector, localVector);
		multiply(localVector, 1-genomes[i]->getWeightFactor());

		// GLOBAL VECTOR CREATION

		std::vector<double> globalVector = genomes[i]->getX();
		std::vector<Genome*> possibilities = genomes;
		possibilities.erase(std::find(possibilities.begin(), possibilities.end(), genomes[i]));

		xr0 = pickRandom(possibilities)->getX();
		xr1 = pickRandom(possibilities)->getX();
		subtract(bestX, genomes[i]->getX(), subtraction);
		multiply(subtraction, alpha);
		add(subtraction, globalVector, globalVector);
		subtract(xr0, xr1, subtraction);
		multiply(subtraction, beta);
		add(subtraction, globalVector, globalVector);
		multiply(globalVector, genomes[i]->getWeightFactor());


		std::vector<double> mutant(D);
		add(localVector, globalVector, mutant);

		std::vector<Genome*> copy = genomes;
		copy.erase(std::find(copy.begin(), copy.end(), genomes[i]));

		Genome* gr0 = pickRandom(copy);
		Genome* gr1 = pickRandom(copy); 

		double newWeight = genomes[i]->getWeightFactor() + (F * (best->getWeightFactor() - genomes[i]->getWeightFactor()))
			+ (F * (gr0->getWeightFactor() - gr1->getWeightFactor()));

		if (newWeight > 0.95)
			newWeight = 0.95;
		else if (newWeight < 0.05)
			newWeight = 0.05;

		Genome * m = new Genome(mutant);
		m->setWeightFactor(newWeight);

		mutants.push_back(m);
	}

	return mutants;
}
