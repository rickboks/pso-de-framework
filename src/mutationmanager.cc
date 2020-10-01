#include "mutationmanager.h"
#include "util.h"
#include <limits>
#include <numeric>
#define LC(X) [](int const D, DEConstraintHandler* const ch){return new X(D,ch);}
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

std::map<std::string, std::function<MutationManager* (int const, DEConstraintHandler*const)>> const mutations ({
		{"R1", LC(Rand1MutationManager)},
		{"T1", LC(TTB1MutationManager)},
		{"T2", LC(TTB2MutationManager)},
		{"P1", LC(TTPB1MutationManager)},
		{"B1", LC(Best1MutationManager)},
		{"B2", LC(Best2MutationManager)},
		{"R2", LC(Rand2MutationManager)},
		{"RD", LC(Rand2DirMutationManager)},
		{"NS", LC(NSDEMutationManager)},
		{"TR", LC(TrigonometricMutationManager)},
		{"O1", LC(TwoOpt1MutationManager)},
		{"O2", LC(TwoOpt2MutationManager)},
		{"PX", LC(ProximityMutationManager)},
});

// Rand/1
Solution* Rand1MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<Solution*> xr = pickRandom(possibilities, 3);
	std::vector<double> mutant = xr[0]->getX();
	std::vector<double> difference(this->D);
	subtract(xr[1]->getX(), xr[2]->getX(), difference);
	scale(difference, Fs[i]);
	add(mutant,difference, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, xr[0], genomes[i]);
	return m;
}

// Target-to-best/1
Solution* TTB1MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);
	std::vector<double> mutant = genomes[i]->getX();
	std::vector<double> difference(this->D);
	std::vector<Solution*> xr = pickRandom(possibilities, 2);

	subtract(best->getX(), genomes[i]->getX(), difference);
	add(difference, xr[0]->getX(), difference);
	subtract(difference, xr[1]->getX(), difference);
	scale(difference,Fs[i]);

	add(mutant, difference, mutant);
	Solution* m = new Solution(mutant);
	deCH->repairDE(m, genomes[i], genomes[i]);
	return m;
}

// Target-to-best/2
Solution* TTB2MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = genomes[i]->getX();
	std::vector<double> difference(this->D);

	std::vector<Solution*> xr = pickRandom(possibilities, 4);

	subtract(best->getX(), genomes[i]->getX(), difference);
	add(difference, xr[0]->getX(), difference);
	subtract(difference, xr[1]->getX(), difference);
	add(difference, xr[2]->getX(), difference);
	subtract(difference, xr[3]->getX(), difference);
	scale(difference,Fs[i]);

	add(mutant, difference, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, genomes[i], genomes[i]);
	return m;
}

Solution* TTPB1MutationManager::mutate(int const i) const{
	Solution* pBest = getPBest(this->genomes);
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = genomes[i]->getX();
	std::vector<double> difference(this->D);

	std::vector<Solution*> xr = pickRandom(possibilities, 4);

	subtract(pBest->getX(), genomes[i]->getX(), difference);
	add(difference, xr[0]->getX(), difference);
	subtract(difference, xr[1]->getX(), difference);
	add(difference, xr[2]->getX(), difference);
	subtract(difference, xr[3]->getX(), difference);
	scale(difference,Fs[i]);

	add(mutant, difference, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, genomes[i], genomes[i]);
	return m;
}

// Best/1
Solution* Best1MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = best->getX();
	std::vector<double> difference(this->D);
	
	std::vector<Solution*> xr = pickRandom(possibilities, 2);
	subtract(xr[0]->getX(), xr[1]->getX(), difference);

	scale(difference, Fs[i]);
	add(mutant, difference, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, best, genomes[i]);
	return m;
}

// Best/2
Solution* Best2MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = best->getX();
	std::vector<double> difference(this->D);

	std::vector<Solution*> xr = pickRandom(possibilities, 4);
	subtract(xr[0]->getX(), xr[1]->getX(), difference);
	add(difference, xr[2]->getX(), difference);
	subtract(difference, xr[3]->getX(), difference);
	scale(difference, Fs[i]);
	add(mutant, difference, mutant);

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
	std::vector<double> difference(this->D);

	subtract(xr[0]->getX(), xr[1]->getX(), difference);
	add(difference, xr[2]->getX(), difference);
	subtract(difference, xr[3]->getX(), difference);
	scale(difference, Fs[i]);

	add(mutant, difference, mutant);
	Solution* m = new Solution(mutant);
	deCH->repairDE(m, xr[4], genomes[i]);
	return m;
}

// Rand/2/dir
Solution* Rand2DirMutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<Solution*> xr = pickRandom(possibilities, 3);
	sortOnFitness(xr);

	std::vector<double> mutant = xr[0]->getX();

	std::vector<double> difference = xr[0]->getX();
	add(difference, xr[0]->getX(), difference);
	subtract(difference, xr[1]->getX(), difference);
	subtract(difference, xr[2]->getX(), difference);
	scale(difference, Fs[i]/2.);

	add(mutant, difference, mutant);

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

	std::vector<double> difference(this->D);
	subtract(xr[1]->getX(), xr[2]->getX(), difference);

	double randomVar;
	if (rng.randDouble(0,1) < 0.5)
		randomVar = rng.normalDistribution(0.5,0.5);
	else 
		randomVar = rng.cauchyDistribution(0,1);

	scale(difference, randomVar);

	add(mutant, difference, mutant);

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

	double const pPrime = std::abs(xr[0]->getFitness()) + std::abs(xr[1]->getFitness()) 
					+ std::abs(xr[2]->getFitness());

	double const p0 = xr[0]->getFitness() / pPrime;
	double const p1 = xr[1]->getFitness() / pPrime;
	double const p2 = xr[2]->getFitness() / pPrime;

	std::vector<double> mutant;

	add(xr[0]->getX(), xr[1]->getX(), mutant);
	add(mutant, xr[2]->getX(), mutant);
	scale(mutant, 1./3.);

	Solution base = Solution(mutant); // only used for correction strategies

	std::vector<double> temp(this->D);
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
	std::vector<double> difference(this->D);
	subtract(xr[1]->getX(), xr[2]->getX(), difference);
	scale(difference, Fs[i]);
	add(mutant,difference, mutant);

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
	std::vector<double> difference(this->D);
	subtract(xr[1]->getX(), xr[2]->getX(), difference);
	scale(difference, Fs[i]);
	add(mutant,difference, mutant);

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
	std::vector<double> difference(this->D);

	subtract(xr[1]->getX(), xr[2]->getX(), difference);
	add(difference, xr[3]->getX(), difference);
	subtract(difference, xr[4]->getX(), difference);
	scale(difference, Fs[i]);

	add(mutant, difference, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, xr[0], genomes[i]);
	return m;
}

// Proximity-based Rand/1
Solution* ProximityMutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> dist(possibilities.size());
	std::vector<double> prob(possibilities.size());

	double totalDist=0;
	for (unsigned int j = 0; j < possibilities.size(); j++){
		dist[j] = distance(genomes[i], possibilities[j]);
		totalDist+=dist[j];
	}

	for (unsigned int j = 0; j < possibilities.size(); j++)
		prob[j] = 1-(dist[j]/std::max(totalDist, std::numeric_limits<double>::min())); //totalDist could be 0

	std::vector<Solution*> xr = rouletteSelect(possibilities, prob, 3);

	std::vector<double> mutant = xr[0]->getX();
	std::vector<double> difference(this->D);
	subtract(xr[1]->getX(), xr[2]->getX(), difference);
	scale(difference, Fs[i]);
	add(mutant,difference, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, xr[0], genomes[i]);
	return m;
}
