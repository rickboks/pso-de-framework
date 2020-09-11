#include "mutationmanager.h"

MutationManager::MutationManager(int const D, ConstraintHandler* deCH)
			:D(D), deCH(deCH){}

MutationManager::~MutationManager(){}

MutationManager* MutationManager::createMutationManager(MutationType const mutationType, int const D, ConstraintHandler* deCH){
	switch(mutationType){
		case RAND_1:
			return new Rand1MutationManager(D, deCH);
		case TTB_1:
			return new TTB1MutationManager(D, deCH);
		case TTPB_1:
			return new TTPB1MutationManager(D, deCH);
		case BEST_1:
			return new Best1MutationManager(D, deCH);
		case BEST_2:
			return new Best2MutationManager(D, deCH);
		case RAND_2:
			return new Rand2MutationManager(D, deCH);
		case RAND_2_DIR:
			return new Rand2DirMutationManager(D, deCH);
		case NSDE:
			return new NSDEMutationManager(D, deCH);
		case TRIGONOMETRIC:
			return new TrigonometricMutationManager(D, deCH);
		case TO1:
			return new TwoOpt1MutationManager(D, deCH);
		case TO2:
			return new TwoOpt2MutationManager(D, deCH);
		// case TOPOLOGY:
		// 	return new TopologyMutationManager(genomes,F,D);
		default:
			throw std::invalid_argument("Error: Invalid DE mutation type");
	}	
}


Rand1MutationManager::Rand1MutationManager(int const D, ConstraintHandler* deCH)
	: MutationManager(D, deCH){}

std::vector<Particle*> Rand1MutationManager::mutate(std::vector<Particle*> const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<Particle*> xr = pickRandom(possibilities, 3);
		std::vector<double> x = xr[0]->getX();
		std::vector<double> subtraction(this->D);
		subtract(xr[1]->getX(), xr[2]->getX(), subtraction);
		scale(subtraction, Fs[i]);
		add(x,subtraction, x);

		Particle* mutant = new Particle(x);
		deCH->repair(mutant, xr[0], genomes[i]);
		mutants.push_back(mutant);		
	}
	return mutants;
}

TTB1MutationManager::TTB1MutationManager(int const D, ConstraintHandler* deCH)
	: MutationManager(D, deCH){}

std::vector<Particle*> TTB1MutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<double> best = getBest(genomes)->getX(); 
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<double> mutant = genomes[i]->getX();
		std::vector<double> subtraction(this->D);
		subtract(best, genomes[i]->getX(), subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);

		std::vector<Particle*> xr = pickRandom(possibilities, 2);
		subtract(xr[0]->getX(), xr[1]->getX(), subtraction);
		scale(subtraction,Fs[i]);
		add(mutant, subtraction, mutant);
		Particle* m = new Particle(mutant);

		deCH->repair(m, genomes[i], genomes[i]);
		mutants.push_back(m);
	}
	
	return mutants;
}

TTPB1MutationManager::TTPB1MutationManager(int const D, ConstraintHandler* deCH)
	:MutationManager(D, deCH), p(0.1){}


std::vector<Particle*> TTPB1MutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<double> best = getPBest(genomes, p)->getX(); 
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<double> mutant = genomes[i]->getX();
		std::vector<double> subtraction(this->D);
		subtract(best, genomes[i]->getX(), subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);

		std::vector<Particle*> xr = pickRandom(possibilities, 2);
		subtract(xr[0]->getX(), xr[1]->getX(), subtraction);
		scale(subtraction,Fs[i]);
		add(mutant, subtraction, mutant);

		Particle* m = new Particle(mutant);
		deCH->repair(m, genomes[i], genomes[i]);
		mutants.push_back(m);
	}
	
	return mutants;
}

Best1MutationManager::Best1MutationManager(int const D, ConstraintHandler* deCH):MutationManager(D, deCH){}
std::vector<Particle*> Best1MutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());
	Particle* best = getBest(genomes);

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<double> mutant = best->getX();
		std::vector<double> subtraction(this->D);
		
		std::vector<Particle*> xr = pickRandom(possibilities, 2);

		subtract(xr[0]->getX(), xr[1]->getX(), subtraction);
		scale(subtraction,Fs[i]);
		add(mutant, subtraction, mutant);

		Particle* m = new Particle(mutant);
		deCH->repair(m, best, genomes[i]);
		mutants.push_back(m);
	}

	return mutants;
}

Best2MutationManager::Best2MutationManager(int const D, ConstraintHandler* deCH):MutationManager(D, deCH){}

std::vector<Particle*> Best2MutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	Particle* best = getBest(genomes);
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<double> mutant = best->getX();
		std::vector<double> subtraction(this->D);

		std::vector<Particle*> xr = pickRandom(possibilities, 4);
		subtract(xr[0]->getX(), xr[1]->getX(), subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);
		subtract(xr[2]->getX(), xr[3]->getX(), subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);

		Particle* m = new Particle(mutant);
		deCH->repair(m, best, genomes[i]);
		mutants.push_back(m);
	}

	return mutants;
}

Rand2MutationManager::Rand2MutationManager(int const D, ConstraintHandler* deCH):MutationManager(D, deCH){}

std::vector<Particle*> Rand2MutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<Particle*> xr = pickRandom(possibilities, 5);
		std::vector<double> mutant = xr[4]->getX();
		std::vector<double> subtraction(this->D);
		subtract(xr[0]->getX(), xr[1]->getX(), subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);
		subtract(xr[2]->getX(), xr[3]->getX(), subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);
		Particle* m = new Particle(mutant);

		deCH->repair(m, xr[4], genomes[i]);
		mutants.push_back(m);
	}
	
	return mutants;
}

Rand2DirMutationManager::Rand2DirMutationManager(int const D, ConstraintHandler* deCH):MutationManager(D, deCH){}

std::vector<Particle*> Rand2DirMutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<Particle*> xr = pickRandom(possibilities, 4);

		if (xr[1]->getFitness() < xr[0]->getFitness()){
			std::swap(xr[0], xr[1]);
		}

		if (xr[3]->getFitness() < xr[2]->getFitness()){
			std::swap(xr[3], xr[2]);
		}

		std::vector<double> mutant = xr[0]->getX();
		std::vector<double> subtraction(this->D);
		subtract(xr[0]->getX(), xr[1]->getX(), subtraction);
		add(subtraction, xr[2]->getX(), subtraction);
		subtract(subtraction, xr[3]->getX(), subtraction);
		scale(subtraction, Fs[i]*0.5);
		add(mutant, subtraction, mutant);

		Particle* m = new Particle(mutant);
		deCH->repair(m, xr[0], genomes[i]);
		mutants.push_back(m);
	}
	
	return mutants;
}

NSDEMutationManager::NSDEMutationManager(int const D, ConstraintHandler* deCH):MutationManager(D, deCH){}

std::vector<Particle*> NSDEMutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);


		std::vector<Particle*> xr = pickRandom(possibilities, 3);
		std::vector<double> mutant = xr[0]->getX(); 

		std::vector<double> subtraction(this->D);
		subtract(xr[1]->getX(), xr[2]->getX(), subtraction);

		if (rng.randDouble(0,1) < 0.5){
			double n = rng.normalDistribution(0.5,0.5);
			scale(subtraction, n);
		} else {
			double c = rng.cauchyDistribution(0,1);
			scale(subtraction, c);
		}

		add(mutant, subtraction, mutant);

		Particle* m = new Particle(mutant);
		deCH->repair(m, xr[0], genomes[i]);
		mutants.push_back(m);
	}

	return mutants;
}

TrigonometricMutationManager::TrigonometricMutationManager(int const D, ConstraintHandler* deCH): MutationManager(D, deCH), gamma(0.05){}
std::vector<Particle*> TrigonometricMutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	return rng.randDouble(0,1) <= gamma ? trigonometricMutation(genomes,Fs) : rand1Mutation(genomes,Fs);
}
std::vector<Particle*> TrigonometricMutationManager::trigonometricMutation(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);
		
		std::vector<Particle*> xr = pickRandom(possibilities, 3);

		double pPrime = fabs(xr[0]->getFitness()) + fabs(xr[1]->getFitness()) 
						+ fabs(xr[2]->getFitness());

		double p0 = xr[0]->getFitness() / pPrime;
		double p1 = xr[1]->getFitness() / pPrime;
		double p2 = xr[2]->getFitness() / pPrime;

		std::vector<double> mutant;
		std::vector<double> temp(this->D);

		add(xr[0]->getX(), xr[1]->getX(), temp);
		add(temp, xr[2]->getX(), temp);
		scale(temp, 1.0/3.0);
		mutant = temp;

		Particle base = Particle(mutant);

		subtract(xr[0]->getX(), xr[1]->getX(), temp);
		scale(temp, p1-p0);
		add(temp, mutant, mutant);

		subtract(xr[1]->getX(), xr[2]->getX(), temp);
		scale(temp, p2-p1);
		add(temp, mutant, mutant);

		subtract(xr[2]->getX(), xr[0]->getX(), temp);
		scale(temp, p0-p2);
		add(temp, mutant, mutant);
		
		Particle* m = new Particle(mutant);
		deCH->repair(m, &base, genomes[i]);
		mutants.push_back(m);
	}

	return mutants;
}

std::vector<Particle*> TrigonometricMutationManager::rand1Mutation(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);
		
		std::vector<Particle*> xr = pickRandom(possibilities, 3);

		std::vector<double> x = xr[0]->getX();
		std::vector<double> subtraction(this->D);
		subtract(xr[1]->getX(), xr[2]->getX(), subtraction);
		scale(subtraction, Fs[i]);
		add(x,subtraction, x);

		Particle* m = new Particle(x);
		deCH->repair(m, xr[0], genomes[i]);
		mutants.push_back(m);
	}
	return mutants;	
}

TwoOpt1MutationManager::TwoOpt1MutationManager(int const D, ConstraintHandler* deCH): MutationManager(D, deCH){}
std::vector<Particle*> TwoOpt1MutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<Particle*> xr = pickRandom(possibilities, 3);

		if (xr[1]->getFitness() < xr[0]->getFitness()){
			std::swap(xr[0], xr[1]);
		}

		std::vector<double> x = xr[0]->getX();
		std::vector<double> subtraction(this->D);
		subtract(xr[1]->getX(), xr[2]->getX(), subtraction);
		scale(subtraction, Fs[i]);
		add(x,subtraction, x);

		Particle* mutant = new Particle(x);
		deCH->repair(mutant, xr[0], genomes[i]);
		mutants.push_back(mutant);		
	}
	return mutants;
}

TwoOpt2MutationManager::TwoOpt2MutationManager(int const D, ConstraintHandler* deCH): MutationManager(D, deCH){}
std::vector<Particle*> TwoOpt2MutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<Particle*> xr = pickRandom(possibilities, 5);

		if (xr[1]->getFitness() < xr[0]->getFitness()){
			std::swap(xr[0], xr[1]);
		}

		std::vector<double> mutant = xr[0]->getX();
		std::vector<double> subtraction(this->D);
		subtract(xr[1]->getX(), xr[2]->getX(), subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);
		subtract(xr[3]->getX(), xr[4]->getX(), subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);

		Particle* m = new Particle(mutant);
		deCH->repair(m, xr[0], genomes[i]);
		mutants.push_back(m);
	}
	
	return mutants;
}
