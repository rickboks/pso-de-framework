#include "selectionmanager.h"
#include "vectoroperations.h"
#include "deadaptationmanager.h"
#include <tuple>
//BASE
SelectionManager::SelectionManager(int const D, DEAdaptationManager* dam)
: D(D), dam(dam){

}

void SelectionManager::selection(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2){
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
: SelectionManager(D, dam){

}

void Pairwise2SelectionManager::select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2){
	for (unsigned int i = 0 ; i < particles.size(); i++){
		//PSO wins
		if (p0[i]->getFitness() < p2[i]->getFitness()){
			particles[i]->setPosition(p0[i]->getPosition(), p0[i]->getFitness());
			particles[i]->setVelocity(p0[i]->getVelocity());
		} else {
			//DE
			std::vector<double> newVelocity(D);
			//reverse engineer velocity
			subtract(p2[i]->getPosition(), particles[i]->getPosition(), newVelocity);
			particles[i]->setVelocity(newVelocity);
			particles[i]->setPosition(p2[i]->getPosition(), p2[i]->getFitness());							
		}
	}
}


//Pairwise 3
Pairwise3SelectionManager::Pairwise3SelectionManager(int const D, DEAdaptationManager* dam)
: SelectionManager(D, dam){

}

void Pairwise3SelectionManager::select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2){

	for (unsigned int i = 0 ; i < particles.size(); i++){
		double Fparticle = particles[i]->getFitness();
		double Fp0 = p0[i]->getFitness();
		double Fp2 = p2[i]->getFitness();

		if (Fparticle < Fp0 && Fparticle < Fp2){
			//Do nothing
		} else if (Fp0 < Fparticle && Fp0 < Fp2){
			particles[i]->setPosition(p0[i]->getPosition(), p0[i]->getFitness());
			particles[i]->setVelocity(p0[i]->getVelocity());
		} else {
			std::vector<double> newVelocity(D);
			//reverse engineer velocity
			subtract(p2[i]->getPosition(), particles[i]->getPosition(), newVelocity);
			particles[i]->setVelocity(newVelocity);
			particles[i]->setPosition(p2[i]->getPosition(), p2[i]->getFitness());
		}
	}
}

//Union 2
Union2SelectionManager::Union2SelectionManager(int const D, DEAdaptationManager* dam)
: SelectionManager(D, dam){

}

bool sortbyFirstElement(const std::tuple<Particle*, int, int>& a,  
               const std::tuple<Particle*, int, int>& b) { 
    return (std::get<0>(a)->getFitness() < std::get<0>(b)->getFitness()); 
} 
  
void Union2SelectionManager::select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2){

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
			particles[index]->setPosition(p->getPosition(), p->getFitness());
			particles[index]->setVelocity(p->getVelocity());
		} else {
			std::vector<double> newVelocity(D);			//reverse engineer velocity
			subtract(p->getPosition(), particles[index]->getPosition(), newVelocity);
			particles[index]->setVelocity(newVelocity);
			particles[index]->setPosition(p->getPosition(), p->getFitness());
		}
	}
}


//Union 3
Union3SelectionManager::Union3SelectionManager(int const D, DEAdaptationManager* dam)
: SelectionManager(D, dam){

}

void Union3SelectionManager::select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2){
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
			particles[index]->setPosition(p->getPosition(), p->getFitness());
			particles[index]->setVelocity(p->getVelocity());
		} else if (type == 1){
			std::vector<double> newVelocity(D);			//reverse engineer velocity
			subtract(p->getPosition(), particles[index]->getPosition(), newVelocity);
			particles[index]->setVelocity(newVelocity);
			particles[index]->setPosition(p->getPosition(), p->getFitness());
		} else {
			// do nothing: original is best
		}
	}

}

SelectionManager* SelectionManagerFactory::createSelectionManager(SelectionType type, int const D, DEAdaptationManager* dam){
	switch(type){
		case P2:
			return new Pairwise2SelectionManager(D, dam);
		case P3:
			return new Pairwise3SelectionManager(D, dam);
		case U2:
			return new Union2SelectionManager(D, dam);
		case U3:
			return new Union3SelectionManager(D, dam);
		default:
			throw std::invalid_argument("Error: Invalid DE mutation type");
	}
}
