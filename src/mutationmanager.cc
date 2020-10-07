#include "mutationmanager.h"
#include "util.h"
#include <limits>
#include <numeric>

#define LC(X) [](int const D, DEConstraintHandler* const ch){return new X(D,ch);}

std::vector<Solution*> MutationManager::mutate(std::vector<Solution*>const& genomes, std::vector<double>const& Fs){
	this->genomes = genomes;
	this->Fs = Fs;

	preMutation(); // Some mutation managers use this to prepare some stuff

	std::vector<Solution*> mutants(genomes.size());

	for (unsigned int i = 0; i < genomes.size(); i++){
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
		{"RA", LC(RankingMutationManager)},
});

// Rand/1
Solution* Rand1MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<Solution*> xr = pickRandom(possibilities, 3);
	std::vector<double> mutant = xr[0]->getX();
	std::vector<double> difference(D);
	subtract(xr[1]->getX(), xr[2]->getX(), difference);
	scale(difference, Fs[i]);
	add(mutant,difference, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, xr[0], genomes[i]);
	return m;
}

// Target-to-best/1
void TTB1MutationManager::preMutation(){
	best = getBest(genomes);
}

Solution* TTB1MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = genomes[i]->getX();
	std::vector<double> difference(D);

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
void TTB2MutationManager::preMutation(){
	best = getBest(genomes);
}

Solution* TTB2MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = genomes[i]->getX();
	std::vector<double> difference(D);

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

// Target-to-pbest/1
Solution* TTPB1MutationManager::mutate(int const i) const{
	Solution* pBest = getPBest(genomes); // pBest is sampled for each mutation

	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = genomes[i]->getX();
	std::vector<double> difference(D);

	std::vector<Solution*> xr = pickRandom(possibilities, 2);

	subtract(pBest->getX(), genomes[i]->getX(), difference);
	add(difference, xr[0]->getX(), difference);
	subtract(difference, xr[1]->getX(), difference);
	scale(difference,Fs[i]);

	add(mutant, difference, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, genomes[i], genomes[i]);
	return m;
}

// Best/1
void Best1MutationManager::preMutation(){
	best = getBest(genomes);
}

Solution* Best1MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = best->getX();
	std::vector<double> difference(D);
	
	std::vector<Solution*> xr = pickRandom(possibilities, 2);
	subtract(xr[0]->getX(), xr[1]->getX(), difference);

	scale(difference, Fs[i]);
	add(mutant, difference, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, best, genomes[i]);
	return m;
}

// Best/2
void Best2MutationManager::preMutation(){
	best = getBest(genomes);
}

Solution* Best2MutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = best->getX();
	std::vector<double> difference(D);

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
	std::vector<double> difference(D);

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

	std::vector<Solution*> xr = pickRandom(possibilities, 4);

	if (xr[1]->getFitness() < xr[0]->getFitness())
		std::swap(xr[0], xr[1]);

	if (xr[3]->getFitness() < xr[2]->getFitness())
		std::swap(xr[2], xr[3]);

	std::vector<double> mutant = xr[0]->getX();

	std::vector<double> difference(D);
	subtract(xr[0]->getX(), xr[1]->getX(), difference);
	add(difference, xr[2]->getX(), difference);
	subtract(difference, xr[3]->getX(), difference);
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

	std::vector<double> difference(D);
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

	double const p0 = std::abs(xr[0]->getFitness()) / pPrime;
	double const p1 = std::abs(xr[1]->getFitness()) / pPrime;
	double const p2 = std::abs(xr[2]->getFitness()) / pPrime;

	std::vector<double> mutant(D);

	add(xr[0]->getX(), xr[1]->getX(), mutant);
	add(mutant, xr[2]->getX(), mutant);
	scale(mutant, 1./3.);

	Solution base = Solution(mutant); // only used for correction strategies

	std::vector<double> temp(D);
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
	std::vector<double> difference(D);
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
	std::vector<double> difference(D);
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
	std::vector<double> difference(D);

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
// TODO: currently, this assumes hyperbox constraints by using 
// eucledian distance.
// TODO: not implemented correctly I think. Contacted an author.
void ProximityMutationManager::preMutation(){
	int const size = genomes.size();
	//Initialize the matrices
	if (Rp.empty()){
		Rp.resize(size, std::vector<double>(size));
		Rd.resize(size, std::vector<double>(size));
	}

	// Fill distance matrix
	std::vector<double> rowTotals(size, 0.);
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){
			if (i != j){
				double const dist = std::max(distance(genomes[i], genomes[j]), 1.0e-12);
				Rd[i][j] = dist;
				Rd[j][i] = dist;
				rowTotals[i] += dist;
			} else {
				Rd[i][j] = 0.;
			}
		}
	}
	
	// Fill probability matrix
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){
			if (i != j){
				double const prob = 1. / (Rd[i][j] / rowTotals[i]);
				Rp[i][j] = prob;
				Rp[j][i] = prob;
			} else {
				Rp[i][j] = 0.;
			}
		}
	}
}

Solution* ProximityMutationManager::mutate(int const i) const{
	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> prob = Rp[i];
	prob.erase(prob.begin() + i); // Remove own probability
	std::vector<Solution*> xr = rouletteSelect(possibilities, prob, 3);

	std::vector<double> mutant = xr[0]->getX();
	std::vector<double> difference(D);
	subtract(xr[1]->getX(), xr[2]->getX(), difference);
	scale(difference, Fs[i]);
	add(mutant,difference, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, xr[0], genomes[i]);
	return m;
}

// Ranking based
void RankingMutationManager::preMutation(){
	int const size = genomes.size();
	probability.clear();

	std::vector<Solution*> sorted = genomes;
	sortOnFitness(sorted);

	for (int i = 0; i < size; i++)
		probability[sorted[i]] = double(size - (i+1)) / double(size);
}

Solution* RankingMutationManager::pickRanked(std::vector<Solution*>& possibilities) const{
	int index;
	do {
		index = rng.randInt(0, possibilities.size()-1);
	} while (rng.randDouble(0,1) > probability.at(possibilities[index]));

	Solution* const pick = possibilities[index];
	possibilities.erase(possibilities.begin() + index);
	return pick;
}

Solution* RankingMutationManager::mutate(int const i) const{
	Solution* pBest = getPBest(genomes); // pBest is sampled for each mutation

	std::vector<Solution*> possibilities = genomes;
	possibilities.erase(possibilities.begin() + i);

	std::vector<double> mutant = genomes[i]->getX();
	std::vector<double> difference(D);

	Solution* xr0 = pickRanked(possibilities); // N.B. Ranked instead of Random

	Solution* xr1 = pickRandom(possibilities);

	subtract(pBest->getX(), genomes[i]->getX(), difference);
	add(difference, xr0->getX(), difference);
	subtract(difference, xr1->getX(), difference);
	scale(difference,Fs[i]);

	add(mutant, difference, mutant);

	Solution* m = new Solution(mutant);
	deCH->repairDE(m, genomes[i], genomes[i]);
	return m;
}
