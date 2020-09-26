#include "selectionmanager.h"
#include "util.h"
#include "deadaptationmanager.h"
#include <tuple>
//BASE
SelectionManager::SelectionManager(int const D, DEAdaptationManager* const dam)
: D(D), dam(dam){
}

#define LC(X) [](int const D, DEAdaptationManager* const adap){return new X(D,adap);}
std::map<std::string, std::function<SelectionManager* (int const, DEAdaptationManager *const)>> const selections ({
		{"P2", LC(Pairwise2SelectionManager)},
		{"P3", LC(Pairwise3SelectionManager)},
		{"U2", LC(Union2SelectionManager)},
		{"U3", LC(Union3SelectionManager)},
});

void SelectionManager::selection(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Solution*>const& p2) const{
	//check sucessful indices
	for (unsigned int i = 0; i < particles.size(); i++){
		if (p2[i]->getFitness() < particles[i]->getFitness()){
			dam->successfulIndex(i);
		}
	}
	//start selection
	select(particles, p0, p2);
}

//Pairwise 2
Pairwise2SelectionManager::Pairwise2SelectionManager(int const D, DEAdaptationManager* dam) 
: SelectionManager(D, dam){}

void Pairwise2SelectionManager::select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Solution*>const& p2) const{
	for (unsigned int i = 0 ; i < particles.size(); i++){
		//PSO wins
		if (p0[i]->getFitness() < p2[i]->getFitness()){
			particles[i]->setX(p0[i]->getX(), p0[i]->getFitness());
			particles[i]->setV(p0[i]->getV());
		} else {
			//DE
			particles[i]->setX(p2[i]->getX(), p2[i]->getFitness());
		}
	}
}


//Pairwise 3
Pairwise3SelectionManager::Pairwise3SelectionManager(int const D, DEAdaptationManager* dam)
: SelectionManager(D, dam){}

void Pairwise3SelectionManager::select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Solution*>const& p2) const{

	for (unsigned int i = 0 ; i < particles.size(); i++){
		double Fparticle = particles[i]->getFitness();
		double Fp0 = p0[i]->getFitness();
		double Fp2 = p2[i]->getFitness();

		if (Fparticle < Fp0 && Fparticle < Fp2){
			//Do nothing
		} else if (Fp0 < Fparticle && Fp0 < Fp2){
			particles[i]->setX(p0[i]->getX(), p0[i]->getFitness());
			particles[i]->setV(p0[i]->getV());
		} else {
			particles[i]->setX(p2[i]->getX(), p2[i]->getFitness());
		}
	}
}

//Union 2
Union2SelectionManager::Union2SelectionManager(int const D, DEAdaptationManager* dam)
: SelectionManager(D, dam){}

bool sortbyFirstElement(const std::tuple<Particle*, int, int>& a,  
               const std::tuple<Particle*, int, int>& b) { 
    return (std::get<0>(a)->getFitness() < std::get<0>(b)->getFitness()); 
} 
  
void Union2SelectionManager::select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Solution*>const& p2) const{

	int popSize = particles.size();

	std::vector< std::tuple<Particle*, int, int > >  allSolutions;
	allSolutions.reserve(popSize*2);

	for (int i = 0; i < popSize; i++){
		allSolutions.push_back(std::make_tuple(p0[i], 0, i));
		allSolutions.push_back(std::make_tuple(p2[i], 1, i));
	}

	std::sort(allSolutions.begin(), allSolutions.end(), sortbyFirstElement);
	
	for (int i = 0; i < popSize; i++){
		Particle* p = std::get<0>(allSolutions[i]);
		int index = std::get<2>(allSolutions[i]);		
		if (std::get<1>(allSolutions[i]) == 0){
			particles[index]->setX(p->getX(), p->getFitness());
			particles[index]->setV(p->getV());
		} else {
			particles[index]->setX(p->getX(), p->getFitness());
		}
	}
}

//Union 3
Union3SelectionManager::Union3SelectionManager(int const D, DEAdaptationManager* dam)
: SelectionManager(D, dam){}

void Union3SelectionManager::select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Solution*>const& p2) const{
	int popSize = particles.size();

	std::vector< std::tuple<Particle*, int, int > >  allSolutions;
	allSolutions.reserve(popSize*2);

	for (int i = 0; i < popSize; i++){
		allSolutions.push_back(std::make_tuple(p0[i], 0, i));
		allSolutions.push_back(std::make_tuple(p2[i], 1, i));
		allSolutions.push_back(std::make_tuple(particles[i], 2, i));
	}

	std::sort(allSolutions.begin(), allSolutions.end(), sortbyFirstElement);

	for (int i = 0; i < popSize; i++){
		Particle* p = std::get<0>(allSolutions[i]);
		int type = std::get<1>(allSolutions[i]);
		int index = std::get<2>(allSolutions[i]);		

		if (type == 0){
			particles[index]->setX(p->getX(), p->getFitness(), false);
			particles[index]->setV(p->getV());
		} else if (type == 1){
			particles[index]->setX(p->getX(), p->getFitness(), true);
		} else {
			// do nothing: original is best
		}
	}
}
