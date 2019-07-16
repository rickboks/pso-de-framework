#include "mutationmanager.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <functional>
#include <limits>

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
:genomes(genomes), F(F), D(genomes[0]->getDimension()), popSize(genomes.size()), generator(randDev()){

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
	std::uniform_int_distribution<int>  distr(0, possibilities.size() -1);
	int const r = distr(generator);
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

		std::transform( xr1.begin(), xr1.end(),
		                xr2.begin(), subtraction.begin(), 
		                std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
	       std::bind(std::multiplies<double>(), std::placeholders::_1, F));

		std::transform( x.begin(), x.end(),
		                subtraction.begin(), x.begin(), 
		                std::plus<double>());

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

		std::transform( best.begin(), best.end(),
		    (genomes[i]->getX()).begin(), subtraction.begin(), 
		    std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
	       std::bind(std::multiplies<double>(), std::placeholders::_1, F));

		std::transform( mutant.begin(), mutant.end(),
		    subtraction.begin(), mutant.begin(), 
		    std::plus<double>());


		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();

		std::transform( xr0.begin(), xr0.end(),
		    xr1.begin(), subtraction.begin(), 
		    std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
	       std::bind(std::multiplies<double>(), std::placeholders::_1, F));

		std::transform( mutant.begin(), mutant.end(),
		    subtraction.begin(), mutant.begin(), 
		    std::plus<double>());

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

		std::transform( xr0.begin(), xr0.end(),
		    xr1.begin(), subtraction.begin(), 
		    std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
	       std::bind(std::multiplies<double>(), std::placeholders::_1, F));

		std::transform( mutant.begin(), mutant.end(),
		    subtraction.begin(), mutant.begin(), 
		    std::plus<double>());

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

		std::transform(xr0.begin(), xr0.end(),
		    xr1.begin(), subtraction.begin(), 
		    std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
	       std::bind(std::multiplies<double>(), std::placeholders::_1, F));

		std::transform(mutant.begin(), mutant.end(),
		    subtraction.begin(), mutant.begin(), 
		    std::plus<double>());

		std::transform( xr2.begin(), xr2.end(),
		    xr3.begin(), subtraction.begin(), 
		    std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
	       std::bind(std::multiplies<double>(), std::placeholders::_1, F));

		std::transform( mutant.begin(), mutant.end(),
		    subtraction.begin(), mutant.begin(), 
		    std::plus<double>());

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


		std::transform(xr0.begin(), xr0.end(),
		    xr1.begin(), subtraction.begin(), 
		    std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
	       std::bind(std::multiplies<double>(), std::placeholders::_1, F));

		std::transform(mutant.begin(), mutant.end(),
		    subtraction.begin(), mutant.begin(), 
		    std::plus<double>());


		std::transform( xr2.begin(), xr2.end(),
		    xr3.begin(), subtraction.begin(), 
		    std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
	       std::bind(std::multiplies<double>(), std::placeholders::_1, F));

		std::transform( mutant.begin(), mutant.end(),
		    subtraction.begin(), mutant.begin(),
		    std::plus<double>());

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

		std::transform(xr0.begin(), xr0.end(),
		    xr1.begin(), subtraction.begin(), 
		    std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(),
		    xr2.begin(), subtraction.begin(), 
		    std::plus<double>());

		std::transform(subtraction.begin(), subtraction.end(),
		    xr3.begin(), subtraction.begin(), 
		    std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
	       std::bind(std::multiplies<double>(), std::placeholders::_1, F*0.5));

		std::transform(mutant.begin(), mutant.end(),
		    subtraction.begin(), mutant.begin(), 
		    std::plus<double>());	

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

	std::uniform_real_distribution<double>  distr(0, 1);


	for (int i = 0; i < popSize; i++){
		std::vector<Genome*> possibilities = genomes;

		auto position = std::find(possibilities.begin(), possibilities.end(), genomes[i]);
		if (position != possibilities.end())
		    possibilities.erase(position);

		std::vector<double> mutant = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();
		std::vector<double> xr2 = pickRandom(possibilities)->getX();

		std::vector<double> subtraction(D);

		std::transform( xr1.begin(), xr1.end(),
			                xr2.begin(), subtraction.begin(), 
			                std::minus<double>());

		if (distr(generator) < 0.5){
			std::normal_distribution<double> N(0.5,0.5);
			double n = N(generator);
			std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
		       std::bind(std::multiplies<double>(), std::placeholders::_1, n));
		} else {
			std::cauchy_distribution<double> cauchy(0,1);
			double c = cauchy(generator);
			std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
		       std::bind(std::multiplies<double>(), std::placeholders::_1, c));
		}

		std::transform(mutant.begin(), mutant.end(),
		    subtraction.begin(), mutant.begin(), 
		    std::plus<double>());


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
	std::uniform_real_distribution<double>  distr(0, 1);
	for (Genome* g : genomes){
		g->setWeightFactor(distr(generator));
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
		std::transform(bestNeighbor.begin(), bestNeighbor.end(),
			                genomes[i]->getX().begin(), subtraction.begin(), 
			                std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
		       std::bind(std::multiplies<double>(), std::placeholders::_1, alpha));

		std::transform( subtraction.begin(), subtraction.end(),
			               	localVector.begin(), localVector.begin(), 
			                std::plus<double>());

		std::transform( xr0.begin(), xr0.end(),
			                xr1.begin(), subtraction.begin(), 
			                std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
		       std::bind(std::multiplies<double>(), std::placeholders::_1, beta));

		std::transform( subtraction.begin(), subtraction.end(),
			               	localVector.begin(), localVector.begin(), 
			                std::plus<double>());	

		std::transform(localVector.begin(), localVector.end(), localVector.begin(),
		       std::bind(std::multiplies<double>(), std::placeholders::_1, 1 - genomes[i]->getWeightFactor()));	

		// GLOBAL VECTOR CREATION

		std::vector<double> globalVector = genomes[i]->getX();
		std::vector<Genome*> possibilities = genomes;
		possibilities.erase(std::find(possibilities.begin(), possibilities.end(), genomes[i]));

		xr0 = pickRandom(possibilities)->getX();
		xr1 = pickRandom(possibilities)->getX();

		std::transform( bestX.begin(), bestX.end(),
			                genomes[i]->getX().begin(), subtraction.begin(), 
			                std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
		       std::bind(std::multiplies<double>(), std::placeholders::_1, alpha));

		std::transform( subtraction.begin(), subtraction.end(),
			               	globalVector.begin(), globalVector.begin(), 
			                std::plus<double>());

		std::transform( xr0.begin(), xr0.end(),
			                xr1.begin(), subtraction.begin(), 
			                std::minus<double>());

		std::transform(subtraction.begin(), subtraction.end(), subtraction.begin(),
		       std::bind(std::multiplies<double>(), std::placeholders::_1, beta));

		std::transform( subtraction.begin(), subtraction.end(),
			               	globalVector.begin(), globalVector.begin(), 
			                std::plus<double>());

		std::transform(globalVector.begin(), globalVector.end(), globalVector.begin(),
		       std::bind(std::multiplies<double>(), std::placeholders::_1, genomes[i]->getWeightFactor()));


		std::vector<double> mutant(D);
		std::transform( localVector.begin(), localVector.end(),
			               	globalVector.begin(), mutant.begin(), 
			                std::plus<double>());


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
