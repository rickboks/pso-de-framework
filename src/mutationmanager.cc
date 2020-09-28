#include "mutationmanager.h"
#define LC(X) [](int const D, DEConstraintHandler const * const ch){return new X(D,ch);}
std::vector<Solution*> MutationManager::mutate(std::vector<Solution*>const& genomes, std::vector<double>const& Fs){
	this->genomes = genomes;
	this->Fs = Fs;
	best = getBest(this->genomes);
	std::vector<Solution*> mutants(this->genomes.size());

	for (unsigned int i = 0; i < this->genomes.size(); i++){
		int resamples = 0;
		while (true){
			Solution* m = mutate(i);
			if (!deCH->resample(m, resamples)){
				deCH->repair(m); //generic repair
				mutants[i] = m; 
				break;
			}
			resamples++;
			delete m;
		}
	}
	return mutants;
}
std::map<std::string, std::function<MutationManager* (int const, DEConstraintHandler const*const)>> const mutations ({
		{"R1", LC(Rand1MutationManager)},
		{"T1", LC(TTB1MutationManager)},
		{"P1", LC(TTPB1MutationManager)},
		{"B1", LC(Best1MutationManager)},
		{"B2", LC(Best2MutationManager)},
		{"R2", LC(Rand2MutationManager)},
		{"RD", LC(Rand2DirMutationManager)},
		{"NS", LC(NSDEMutationManager)},
		{"TR", LC(TrigonometricMutationManager)},
		{"O1", LC(TwoOpt1MutationManager)},
		{"O2", LC(TwoOpt2MutationManager)},
});

// Rand/1
Solution* Rand1MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<Solution*> xr = pickRandom(possibilities, 3);
	std::vector<double> mutant = xr[0]->getX();
	std::vector<double> subtraction(this->D);
	subtract(xr[1]->getX(), xr[2]->getX(), subtraction);
	scale(subtraction, Fs[i]);
	add(mutant,subtraction, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, xr[0], genomes[i]);
	return m;
}

// Target-to-best/1
Solution* TTB1MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = genomes[i]->getX();
	std::vector<double> subtraction(this->D);
	subtract(best->getX(), genomes[i]->getX(), subtraction);
	scale(subtraction, Fs[i]);
	add(mutant, subtraction, mutant);

	std::vector<Solution*> xr = pickRandom(possibilities, 2);
	subtract(xr[0]->getX(), xr[1]->getX(), subtraction);
	scale(subtraction,Fs[i]);
	add(mutant, subtraction, mutant);
	Solution* m = new Solution(mutant);
	deCH->repairDE(m, genomes[i], genomes[i]);
	return m;
}

Solution* TTPB1MutationManager::mutate(int const i) const{
	Solution* pBest = getPBest(this->genomes);
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = genomes[i]->getX();
	std::vector<double> subtraction(this->D);
	subtract(pBest->getX(), genomes[i]->getX(), subtraction);
	scale(subtraction, Fs[i]);
	add(mutant, subtraction, mutant);

	std::vector<Solution*> xr = pickRandom(possibilities, 2);
	subtract(xr[0]->getX(), xr[1]->getX(), subtraction);
	scale(subtraction,Fs[i]);
	add(mutant, subtraction, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, genomes[i], genomes[i]);
	return m;
}

// Best/1
Solution* Best1MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = best->getX();
	std::vector<double> subtraction(this->D);
	
	std::vector<Solution*> xr = pickRandom(possibilities, 2);
	subtract(xr[0]->getX(), xr[1]->getX(), subtraction);

	scale(subtraction, Fs[i]);
	add(mutant, subtraction, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, best, genomes[i]);
	return m;
}

// Best/2
Solution* Best2MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = best->getX();
	std::vector<double> subtraction(this->D);

	std::vector<Solution*> xr = pickRandom(possibilities, 4);
	subtract(xr[0]->getX(), xr[1]->getX(), subtraction);
	scale(subtraction, Fs[i]);
	add(mutant, subtraction, mutant);
	subtract(xr[2]->getX(), xr[3]->getX(), subtraction);
	scale(subtraction, Fs[i]);
	add(mutant, subtraction, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, best, genomes[i]);
	return m;
}

// Rand/2
Solution* Rand2MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<Solution*> xr = pickRandom(possibilities, 5);
	std::vector<double> mutant = xr[4]->getX();
	std::vector<double> subtraction(this->D);
	subtract(xr[0]->getX(), xr[1]->getX(), subtraction);
	scale(subtraction, Fs[i]);
	add(mutant, subtraction, mutant);
	subtract(xr[2]->getX(), xr[3]->getX(), subtraction);
	scale(subtraction, Fs[i]);
	add(mutant, subtraction, mutant);
	Solution* m = new Solution(mutant);

	deCH->repairDE(m, xr[4], genomes[i]);
	return m;
}

// Rand/2/dir
Solution* Rand2DirMutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<Solution*> xr = pickRandom(possibilities, 4);

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

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, xr[0], genomes[i]);
	return m;
}

// NSDE
Solution* NSDEMutationManager::mutate(int const i) const {
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<Solution*> xr = pickRandom(possibilities, 3);
	std::vector<double> mutant = xr[0]->getX(); 

	std::vector<double> subtraction(this->D);
	subtract(xr[1]->getX(), xr[2]->getX(), subtraction);

	double randomVar;
	if (rng.randDouble(0,1) < 0.5)
		randomVar = rng.normalDistribution(0.5,0.5);
	else 
		randomVar = rng.cauchyDistribution(0,1);

	scale(subtraction, randomVar);

	add(mutant, subtraction, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, xr[0], genomes[i]);
	return m;
}

// Trigonometric
Solution* TrigonometricMutationManager::mutate(int const i) const{
	return rng.randDouble(0,1) <= gamma ? trigonometricMutation(i) : rand1Mutation(i);
}

Solution* TrigonometricMutationManager::trigonometricMutation(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);
	
	std::vector<Solution*> xr = pickRandom(possibilities, 3);

	double pPrime = std::abs(xr[0]->getFitness()) + std::abs(xr[1]->getFitness()) 
					+ std::abs(xr[2]->getFitness());

	double p0 = xr[0]->getFitness() / pPrime;
	double p1 = xr[1]->getFitness() / pPrime;
	double p2 = xr[2]->getFitness() / pPrime;

	std::vector<double> mutant;
	std::vector<double> temp(this->D);

	add(xr[0]->getX(), xr[1]->getX(), temp);
	add(temp, xr[2]->getX(), temp);
	scale(temp, 1.0/3.0);
	mutant = temp;

	Solution base = Solution(mutant); // only used for correction strategies

	subtract(xr[0]->getX(), xr[1]->getX(), temp);
	scale(temp, p1-p0);
	add(temp, mutant, mutant);

	subtract(xr[1]->getX(), xr[2]->getX(), temp);
	scale(temp, p2-p1);
	add(temp, mutant, mutant);

	subtract(xr[2]->getX(), xr[0]->getX(), temp);
	scale(temp, p0-p2);
	add(temp, mutant, mutant);
	
	Solution* m = new Solution(mutant);
	deCH->repairDE(m, &base, genomes[i]);
	return m;
}

Solution* TrigonometricMutationManager::rand1Mutation(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);
	
	std::vector<Solution*> xr = pickRandom(possibilities, 3);

	std::vector<double> mutant = xr[0]->getX();
	std::vector<double> subtraction(this->D);
	subtract(xr[1]->getX(), xr[2]->getX(), subtraction);
	scale(subtraction, Fs[i]);
	add(mutant,subtraction, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, xr[0], genomes[i]);
	return m;
}

// Two-opt/1
Solution* TwoOpt1MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<Solution*> xr = pickRandom(possibilities, 3);

	if (xr[1]->getFitness() < xr[0]->getFitness())
		std::swap(xr[0], xr[1]);

	std::vector<double> mutant = xr[0]->getX();
	std::vector<double> subtraction(this->D);
	subtract(xr[1]->getX(), xr[2]->getX(), subtraction);
	scale(subtraction, Fs[i]);
	add(mutant,subtraction, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, xr[0], genomes[i]);
	return m;
}

// Two-opt/2
Solution* TwoOpt2MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<Solution*> xr = pickRandom(possibilities, 5);

	if (xr[1]->getFitness() < xr[0]->getFitness())
		std::swap(xr[0], xr[1]);

	std::vector<double> mutant = xr[0]->getX();
	std::vector<double> subtraction(this->D);
	subtract(xr[1]->getX(), xr[2]->getX(), subtraction);
	scale(subtraction, Fs[i]);
	add(mutant, subtraction, mutant);
	subtract(xr[3]->getX(), xr[4]->getX(), subtraction);
	scale(subtraction, Fs[i]);
	add(mutant, subtraction, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, xr[0], genomes[i]);
	return m;
}
