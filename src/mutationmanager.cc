#include "mutationmanager.h"

MutationManager::MutationManager(int const D)
			:D(D){
}

MutationManager::~MutationManager(){}

MutationManager* MutationManager::createMutationManager(MutationType const mutationType, int const D){
	switch(mutationType){
		case RAND_1:
			return new Rand1MutationManager(D);
		case TTB_1:
			return new TTB1MutationManager(D);
		case TTPB_1:
			return new TTPB1MutationManager(D);
		case BEST_1:
			return new Best1MutationManager(D);
		case BEST_2:
			return new Best2MutationManager(D);
		case RAND_2:
			return new Rand2MutationManager(D);
		case RAND_2_DIR:
			return new Rand2DirMutationManager(D);
		case NSDE:
			return new NSDEMutationManager(D);
		case TRIGONOMETRIC:
			return new TrigonometricMutationManager(D);
		case TO1:
			return new TwoOpt1MutationManager(D);
		case TO2:
			return new TwoOpt2MutationManager(D);
		// case TOPOLOGY:
		// 	return new TopologyMutationManager(genomes,F,D);
		default:
			throw std::invalid_argument("Error: Invalid DE mutation type");
	}	
}


Rand1MutationManager::Rand1MutationManager(int const D)
	: MutationManager(D){}


std::vector<Particle*> Rand1MutationManager::mutate(std::vector<Particle*> const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();
		std::vector<double> xr2 = pickRandom(possibilities)->getX();

		std::vector<double> x = xr0;
		std::vector<double> subtraction(this->D);
		subtract(xr1, xr2, subtraction);
		scale(subtraction, Fs[i]);
		add(x,subtraction, x);

		Particle* mutant = new Particle(x);
		mutants.push_back(mutant);		
	}
	return mutants;
}

TTB1MutationManager::TTB1MutationManager(int const D):
	MutationManager(D){}

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

		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();

		subtract(xr0, xr1, subtraction);
		scale(subtraction,Fs[i]);
		add(mutant, subtraction, mutant);

		Particle* m = new Particle(mutant);
		mutants.push_back(m);
	}
	
	return mutants;
}

TTPB1MutationManager::TTPB1MutationManager(int const D)
	:MutationManager(D), p(0.1){}


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

		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();

		subtract(xr0, xr1, subtraction);
		scale(subtraction,Fs[i]);
		add(mutant, subtraction, mutant);

		Particle* m = new Particle(mutant);
		mutants.push_back(m);
	}
	
	return mutants;
}

Best1MutationManager::Best1MutationManager(int const D):MutationManager(D){}
std::vector<Particle*> Best1MutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());
	std::vector<double> best = getBest(genomes)->getX();

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<double> mutant = best;
		std::vector<double> subtraction(this->D);

		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();

		subtract(xr0, xr1, subtraction);
		scale(subtraction,Fs[i]);
		add(mutant, subtraction, mutant);

		Particle* m = new Particle(mutant);
		mutants.push_back(m);
	}


	return mutants;
}

Best2MutationManager::Best2MutationManager(int const D):MutationManager(D){}

std::vector<Particle*> Best2MutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<double> best = getBest(genomes)->getX(); 
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<double> mutant = best;
		std::vector<double> subtraction(this->D);

		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();
		std::vector<double> xr2 = pickRandom(possibilities)->getX();
		std::vector<double> xr3 = pickRandom(possibilities)->getX();
		subtract(xr0, xr1, subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);
		subtract(xr2, xr3, subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);

		Particle* m = new Particle(mutant);

		mutants.push_back(m);
	}

	return mutants;
}

Rand2MutationManager::Rand2MutationManager(int const D):MutationManager(D){}

std::vector<Particle*> Rand2MutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();
		std::vector<double> xr2 = pickRandom(possibilities)->getX();
		std::vector<double> xr3 = pickRandom(possibilities)->getX();
		std::vector<double> xr4 = pickRandom(possibilities)->getX();

		std::vector<double> mutant = xr4;
		std::vector<double> subtraction(this->D);
		subtract(xr0, xr1, subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);
		subtract(xr2, xr3, subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);
		Particle* m = new Particle(mutant);

		mutants.push_back(m);
	}
	
	return mutants;
}

Rand2DirMutationManager::Rand2DirMutationManager(int const D):MutationManager(D){}

std::vector<Particle*> Rand2DirMutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		Particle* r0 = pickRandom(possibilities);
		Particle* r1 = pickRandom(possibilities);
		Particle* r2 = pickRandom(possibilities);
		Particle* r3 = pickRandom(possibilities);

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
		std::vector<double> subtraction(this->D);
		subtract(xr0, xr1, subtraction);
		add(subtraction, xr2, subtraction);
		subtract(subtraction, xr3, subtraction);
		scale(subtraction, Fs[i]*0.5);
		add(mutant, subtraction, mutant);

		Particle* m = new Particle(mutant);

		mutants.push_back(m);
	}
	
	return mutants;
}

NSDEMutationManager::NSDEMutationManager(int const D):MutationManager(D){}

std::vector<Particle*> NSDEMutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<double> mutant = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();
		std::vector<double> xr2 = pickRandom(possibilities)->getX();

		std::vector<double> subtraction(this->D);
		subtract(xr1, xr2, subtraction);

		if (rng.randDouble(0,1) < 0.5){
			double n = rng.normalDistribution(0.5,0.5);
			scale(subtraction, n);
		} else {
			double c = rng.cauchyDistribution(0,1);
			scale(subtraction, c);
		}

		add(mutant, subtraction, mutant);

		Particle* m = new Particle(mutant);

		mutants.push_back(m);
	}

	return mutants;
}

TrigonometricMutationManager::TrigonometricMutationManager(int const D): MutationManager(D), gamma(0.05){}

std::vector<Particle*> TrigonometricMutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	return rng.randDouble(0,1) <= gamma ? trigonometricMutation(genomes,Fs) : rand1Mutation(genomes,Fs);
}

std::vector<Particle*> TrigonometricMutationManager::trigonometricMutation(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		Particle* r0 = pickRandom(possibilities);
		Particle* r1 = pickRandom(possibilities);
		Particle* r2 = pickRandom(possibilities);

		double pPrime = fabs(r0->getFitness()) + fabs(r1->getFitness()) 
						+ fabs(r2->getFitness());

		double p0 = r0->getFitness() / pPrime;
		double p1 = r1->getFitness() / pPrime;
		double p2 = r2->getFitness() / pPrime;

		std::vector<double> mutant;
		std::vector<double> temp(this->D);

		add(r0->getX(), r1->getX(), temp);
		add(temp, r2->getX(), temp);
		scale(temp, 1.0/3.0);
		mutant = temp;

		subtract(r0->getX(), r1->getX(), temp);
		scale(temp, p1-p0);
		add(temp, mutant, mutant);

		subtract(r1->getX(), r2->getX(), temp);
		scale(temp, p2-p1);
		add(temp, mutant, mutant);

		subtract(r2->getX(), r0->getX(), temp);
		scale(temp, p0-p2);
		add(temp, mutant, mutant);

		mutants.push_back(new Particle(mutant));
	}


	return mutants;
}

std::vector<Particle*> TrigonometricMutationManager::rand1Mutation(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<double> xr0 = pickRandom(possibilities)->getX();
		std::vector<double> xr1 = pickRandom(possibilities)->getX();
		std::vector<double> xr2 = pickRandom(possibilities)->getX(); 				

		std::vector<double> x = xr0;
		std::vector<double> subtraction(this->D);
		subtract(xr1, xr2, subtraction);
		scale(subtraction, Fs[i]);
		add(x,subtraction, x);

		Particle* mutant = new Particle(x);
		mutants.push_back(mutant);
	}
	return mutants;	
}

TwoOpt1MutationManager::TwoOpt1MutationManager(int const D): MutationManager(D){}

std::vector<Particle*> TwoOpt1MutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		Particle* r0 = pickRandom(possibilities);
		Particle* r1 = pickRandom(possibilities);

		std::vector<double> xr0;
		std::vector<double> xr1;
		std::vector<double> xr2 = pickRandom(possibilities)->getX();

		if (r0->getFitness() < r1->getFitness()){
			xr0 = r0->getX();
			xr1 = r1->getX();
		} else {
			xr1 = r0->getX();
			xr0 = r1->getX();
		}

		std::vector<double> x = xr0;
		std::vector<double> subtraction(this->D);
		subtract(xr1, xr2, subtraction);
		scale(subtraction, Fs[i]);
		add(x,subtraction, x);

		Particle* mutant = new Particle(x);
		mutants.push_back(mutant);		
	}
	return mutants;
}

TwoOpt2MutationManager::TwoOpt2MutationManager(int const D): MutationManager(D){}

std::vector<Particle*> TwoOpt2MutationManager::mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs){
	std::vector<Particle*> mutants;
	mutants.reserve(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
		std::vector<Particle*> possibilities = genomes;
		possibilities.erase(possibilities.begin() + i);

		std::vector<double> xr0;
		std::vector<double> xr1;
		std::vector<double> xr2 = pickRandom(possibilities)->getX();
		std::vector<double> xr3 = pickRandom(possibilities)->getX();
		std::vector<double> xr4 = pickRandom(possibilities)->getX();

		Particle* r0 = pickRandom(possibilities);
		Particle* r1 = pickRandom(possibilities);

		if (r0->getFitness() < r1->getFitness()){
			xr0 = r0->getX();
			xr1 = r1->getX();
		} else {
			xr1 = r0->getX();
			xr0 = r1->getX();
		}

		std::vector<double> mutant = xr4;
		std::vector<double> subtraction(this->D);
		subtract(xr0, xr1, subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);
		subtract(xr2, xr3, subtraction);
		scale(subtraction, Fs[i]);
		add(mutant, subtraction, mutant);
		Particle* m = new Particle(mutant);

		mutants.push_back(m);
	}
	
	return mutants;
}
