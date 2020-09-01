#pragma once
#include <algorithm>
#include <iostream>
#include <functional>
#include <limits>
#include "rng.h"
#include "util.h"
#include "genome.h"

enum MutationType {
	RAND_1,
	BEST_1,
	TTB_1,
	TTPB_1,
	BEST_2,
	RAND_2,
	RAND_2_DIR,
	NSDE,
	TRIGONOMETRIC,
	TO1,
	TO2,
	MUT_END,

	//Unused
	TOPOLOGY
};

template<class T>
class Rand1MutationManager;
template<class T>
class TTB1MutationManager;
template<class T>
class TTPB1MutationManager;
template<class T>
class Best1MutationManager;
template<class T>
class Best2MutationManager;
template<class T>
class Rand2MutationManager;
template<class T>
class Rand2DirMutationManager;
template<class T>
class NSDEMutationManager;
template<class T>
class TrigonometricMutationManager;
template<class T>
class TwoOpt1MutationManager;
template<class T>
class TwoOpt2MutationManager;

template<class T>
class MutationManager {
	protected:
		int const D;
	public:
		
		MutationManager(int const D)
			:D(D){
		}

		virtual ~MutationManager(){

		}

		virtual std::vector<T*> mutate(std::vector<T*>const& genomes, std::vector<double>& Fs) = 0;
		//template<class T>
		static MutationManager<T>* createMutationManager(MutationType const mutationType, int const D){
			switch(mutationType){
				case RAND_1:
					return new Rand1MutationManager<T>(D);
				case TTB_1:
					return new TTB1MutationManager<T>(D);
				case TTPB_1:
					return new TTPB1MutationManager<T>(D);
				case BEST_1:
					return new Best1MutationManager<T>(D);
				case BEST_2:
					return new Best2MutationManager<T>(D);
				case RAND_2:
					return new Rand2MutationManager<T>(D);
				case RAND_2_DIR:
					return new Rand2DirMutationManager<T>(D);
				case NSDE:
					return new NSDEMutationManager<T>(D);
				case TRIGONOMETRIC:
					return new TrigonometricMutationManager<T>(D);
				case TO1:
					return new TwoOpt1MutationManager<T>(D);
				case TO2:
					return new TwoOpt2MutationManager<T>(D);
				// case TOPOLOGY:
				// 	return new TopologyMutationManager(genomes,F,D);
				default:
					throw std::invalid_argument("Error: Invalid DE mutation type");
			}	
		}
};

template<class T>
class Rand1MutationManager : public MutationManager<T> {
	public:
		Rand1MutationManager(int const D):MutationManager<T>(D){}

		std::vector<T*> mutate(std::vector<T*> const& genomes, std::vector<double>& Fs){
			std::vector<T*> mutants;
			mutants.reserve(genomes.size());

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<T*> possibilities = genomes;
				possibilities.erase(possibilities.begin() + i);

				std::vector<double> xr0 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr1 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr2 = pickRandom(possibilities)->getPosition();

				std::vector<double> x = xr0;
				std::vector<double> subtraction(this->D);
				subtract(xr1, xr2, subtraction);
				scale(subtraction, Fs[i]);
				add(x,subtraction, x);

				T* mutant = new T(x);
				mutants.push_back(mutant);		
			}
			return mutants;
		}
};

template<class T>
class TTB1MutationManager : public MutationManager<T> {
	public:
		TTB1MutationManager(int const D):MutationManager<T>(D){}
		std::vector<T*> mutate(std::vector<T*>const& genomes, std::vector<double>& Fs){
			std::vector<double> best = getBest(genomes)->getPosition(); 
			std::vector<T*> mutants;
			mutants.reserve(genomes.size());

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<T*> possibilities = genomes;
				possibilities.erase(possibilities.begin() + i);

				std::vector<double> mutant = genomes[i]->getPosition();
				std::vector<double> subtraction(this->D);
				subtract(best, genomes[i]->getPosition(), subtraction);
				scale(subtraction, Fs[i]);
				add(mutant, subtraction, mutant);

				std::vector<double> xr0 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr1 = pickRandom(possibilities)->getPosition();

				subtract(xr0, xr1, subtraction);
				scale(subtraction,Fs[i]);
				add(mutant, subtraction, mutant);

				T* m = new T(mutant);
				mutants.push_back(m);
			}
			
			return mutants;
		}
};

template<class T>
class TTPB1MutationManager : public MutationManager<T> {
	private:
		double const p;

	public:
		TTPB1MutationManager(int const D)
			:MutationManager<T>(D), p(0.1){}

		std::vector<T*> mutate(std::vector<T*>const& genomes, std::vector<double>& Fs){
			std::vector<double> best = getPBest(genomes, p)->getPosition(); 
			std::vector<T*> mutants;
			mutants.reserve(genomes.size());

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<T*> possibilities = genomes;
				possibilities.erase(possibilities.begin() + i);

				std::vector<double> mutant = genomes[i]->getPosition();
				std::vector<double> subtraction(this->D);
				subtract(best, genomes[i]->getPosition(), subtraction);
				scale(subtraction, Fs[i]);
				add(mutant, subtraction, mutant);

				std::vector<double> xr0 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr1 = pickRandom(possibilities)->getPosition();

				subtract(xr0, xr1, subtraction);
				scale(subtraction,Fs[i]);
				add(mutant, subtraction, mutant);

				T* m = new T(mutant);
				mutants.push_back(m);
			}
			
			return mutants;
		}
};


template<class T>
class Best1MutationManager: public MutationManager<T> {
	private:

	public:
		Best1MutationManager(int const D):MutationManager<T>(D){}
		std::vector<T*> mutate(std::vector<T*>const& genomes, std::vector<double>& Fs){
			std::vector<T*> mutants;
			mutants.reserve(genomes.size());
			std::vector<double> best = getBest(genomes)->getPosition();

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<T*> possibilities = genomes;
				possibilities.erase(possibilities.begin() + i);

				std::vector<double> mutant = best;
				std::vector<double> subtraction(this->D);

				std::vector<double> xr0 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr1 = pickRandom(possibilities)->getPosition();

				subtract(xr0, xr1, subtraction);
				scale(subtraction,Fs[i]);
				add(mutant, subtraction, mutant);

				T* m = new T(mutant);
				mutants.push_back(m);
			}


			return mutants;
		}
};

template<class T>
class Best2MutationManager: public MutationManager<T> {
	private:

	public:
		Best2MutationManager(int const D):MutationManager<T>(D){}
		std::vector<T*> mutate(std::vector<T*>const& genomes, std::vector<double>& Fs){
			std::vector<double> best = getBest(genomes)->getPosition(); 
			std::vector<T*> mutants;
			mutants.reserve(genomes.size());

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<T*> possibilities = genomes;
				possibilities.erase(possibilities.begin() + i);

				std::vector<double> mutant = best;
				std::vector<double> subtraction(this->D);

				std::vector<double> xr0 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr1 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr2 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr3 = pickRandom(possibilities)->getPosition();
				subtract(xr0, xr1, subtraction);
				scale(subtraction, Fs[i]);
				add(mutant, subtraction, mutant);
				subtract(xr2, xr3, subtraction);
				scale(subtraction, Fs[i]);
				add(mutant, subtraction, mutant);

				T* m = new T(mutant);

				mutants.push_back(m);
			}

			return mutants;
		}
};

template<class T>
class Rand2MutationManager: public MutationManager<T> {
	private:

	public:
		Rand2MutationManager(int const D):MutationManager<T>(D){}
		std::vector<T*> mutate(std::vector<T*>const& genomes, std::vector<double>& Fs){
			std::vector<T*> mutants;
			mutants.reserve(genomes.size());

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<T*> possibilities = genomes;
				possibilities.erase(possibilities.begin() + i);

				std::vector<double> xr0 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr1 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr2 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr3 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr4 = pickRandom(possibilities)->getPosition();

				std::vector<double> mutant = xr4;
				std::vector<double> subtraction(this->D);
				subtract(xr0, xr1, subtraction);
				scale(subtraction, Fs[i]);
				add(mutant, subtraction, mutant);
				subtract(xr2, xr3, subtraction);
				scale(subtraction, Fs[i]);
				add(mutant, subtraction, mutant);
				T* m = new T(mutant);

				mutants.push_back(m);
			}
			
			return mutants;
		}
};

template<class T>
class Rand2DirMutationManager : public MutationManager<T> {
	private:

	public:
		Rand2DirMutationManager(int const D):MutationManager<T>(D){}
		std::vector<T*> mutate(std::vector<T*>const& genomes, std::vector<double>& Fs){
			std::vector<T*> mutants;
			mutants.reserve(genomes.size());

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<T*> possibilities = genomes;
				possibilities.erase(possibilities.begin() + i);

				T* r0 = pickRandom(possibilities);
				T* r1 = pickRandom(possibilities);
				T* r2 = pickRandom(possibilities);
				T* r3 = pickRandom(possibilities);

				std::vector<double> xr0, xr1, xr2, xr3;

				if (r0->getFitness() < r1->getFitness()){
					xr0 = r0->getPosition();
					xr1 = r1->getPosition();
				} else {
					xr0 = r1->getPosition();
					xr1 = r0->getPosition();
				}

				if (r2->getFitness() < r3->getFitness()){
					xr2 = r2->getPosition();
					xr3 = r3->getPosition();
				} else {
					xr2 = r3->getPosition();
					xr3 = r2->getPosition();
				}

				std::vector<double> mutant = xr0;
				std::vector<double> subtraction(this->D);
				subtract(xr0, xr1, subtraction);
				add(subtraction, xr2, subtraction);
				subtract(subtraction, xr3, subtraction);
				scale(subtraction, Fs[i]*0.5);
				add(mutant, subtraction, mutant);

				T* m = new T(mutant);

				mutants.push_back(m);
			}
			
			return mutants;
		}
};

template<class T>
class NSDEMutationManager : public MutationManager<T> {
	private:

	public:
		NSDEMutationManager(int const D):MutationManager<T>(D){}
		std::vector<T*> mutate(std::vector<T*>const& genomes, std::vector<double>& Fs){
			std::vector<T*> mutants;
			mutants.reserve(genomes.size());

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<T*> possibilities = genomes;
				possibilities.erase(possibilities.begin() + i);

				std::vector<double> mutant = pickRandom(possibilities)->getPosition();
				std::vector<double> xr1 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr2 = pickRandom(possibilities)->getPosition();

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

				T* m = new T(mutant);

				mutants.push_back(m);
			}

			return mutants;
		}
};

// class TopologyMutationManager : public MutationManager {
// 	private:
// 		int const radius;
// 		double const alpha;
// 		double const beta;
// 		std::vector<T*> getNeighbors(int const i) const;
// 	public:
// 		TopologyMutationManager(std::vector<double>& Fs, int const D);
// 		std::vector<T*> mutate(std::vector<T*> genomes);
// };

template<class T>
class TrigonometricMutationManager : public MutationManager<T> {
	private:
		double const gamma;
		std::vector<T*> trigonometricMutation(std::vector<T*>const& genomes, std::vector<double>& Fs){
			std::vector<T*> mutants;
			mutants.reserve(genomes.size());

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<T*> possibilities = genomes;
				possibilities.erase(possibilities.begin() + i);

				T* r0 = pickRandom(possibilities);
				T* r1 = pickRandom(possibilities);
				T* r2 = pickRandom(possibilities);

				double pPrime = fabs(r0->getFitness()) + fabs(r1->getFitness()) 
								+ fabs(r2->getFitness());

				double p0 = r0->getFitness() / pPrime;
				double p1 = r1->getFitness() / pPrime;
				double p2 = r2->getFitness() / pPrime;

				std::vector<double> mutant;
				std::vector<double> temp(this->D);

				add(r0->getPosition(), r1->getPosition(), temp);
				add(temp, r2->getPosition(), temp);
				scale(temp, 1.0/3.0);
				mutant = temp;

				subtract(r0->getPosition(), r1->getPosition(), temp);
				scale(temp, p1-p0);
				add(temp, mutant, mutant);

				subtract(r1->getPosition(), r2->getPosition(), temp);
				scale(temp, p2-p1);
				add(temp, mutant, mutant);

				subtract(r2->getPosition(), r0->getPosition(), temp);
				scale(temp, p0-p2);
				add(temp, mutant, mutant);

				mutants.push_back(new T(mutant));
			}


			return mutants;
		}

		std::vector<T*> rand1Mutation(std::vector<T*>const& genomes, std::vector<double>& Fs){
			std::vector<T*> mutants;
			mutants.reserve(genomes.size());

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<T*> possibilities = genomes;
				possibilities.erase(possibilities.begin() + i);

				std::vector<double> xr0 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr1 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr2 = pickRandom(possibilities)->getPosition(); 				

				std::vector<double> x = xr0;
				std::vector<double> subtraction(this->D);
				subtract(xr1, xr2, subtraction);
				scale(subtraction, Fs[i]);
				add(x,subtraction, x);

				T* mutant = new T(x);
				mutants.push_back(mutant);
			}
			return mutants;	
		}

	public:
		TrigonometricMutationManager(int const D): MutationManager<T>(D), gamma(0.05){}
		std::vector<T*> mutate(std::vector<T*>const& genomes, std::vector<double>& Fs){
			return rng.randDouble(0,1) <= gamma ? trigonometricMutation(genomes,Fs) : rand1Mutation(genomes,Fs);
		}
};

template<class T>
class TwoOpt1MutationManager : public MutationManager<T> {
	private:
	public:
		TwoOpt1MutationManager(int const D): MutationManager<T>(D){}
		std::vector<T*> mutate(std::vector<T*>const& genomes, std::vector<double>& Fs){
			std::vector<T*> mutants;
			mutants.reserve(genomes.size());

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<T*> possibilities = genomes;
				possibilities.erase(possibilities.begin() + i);

				T* r0 = pickRandom(possibilities);
				T* r1 = pickRandom(possibilities);

				std::vector<double> xr0;
				std::vector<double> xr1;
				std::vector<double> xr2 = pickRandom(possibilities)->getPosition();

				if (r0->getFitness() < r1->getFitness()){
					xr0 = r0->getPosition();
					xr1 = r1->getPosition();
				} else {
					xr1 = r0->getPosition();
					xr0 = r1->getPosition();
				}

				std::vector<double> x = xr0;
				std::vector<double> subtraction(this->D);
				subtract(xr1, xr2, subtraction);
				scale(subtraction, Fs[i]);
				add(x,subtraction, x);

				T* mutant = new T(x);
				mutants.push_back(mutant);		
			}
			return mutants;
		}
};

template<class T>
class TwoOpt2MutationManager : public MutationManager<T> {
	private:
	public:
		TwoOpt2MutationManager(int const D): MutationManager<T>(D){}
		std::vector<T*> mutate(std::vector<T*>const& genomes, std::vector<double>& Fs){
			std::vector<T*> mutants;
			mutants.reserve(genomes.size());

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<T*> possibilities = genomes;
				possibilities.erase(possibilities.begin() + i);

				std::vector<double> xr0;
				std::vector<double> xr1;
				std::vector<double> xr2 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr3 = pickRandom(possibilities)->getPosition();
				std::vector<double> xr4 = pickRandom(possibilities)->getPosition();

				T* r0 = pickRandom(possibilities);
				T* r1 = pickRandom(possibilities);

				if (r0->getFitness() < r1->getFitness()){
					xr0 = r0->getPosition();
					xr1 = r1->getPosition();
				} else {
					xr1 = r0->getPosition();
					xr0 = r1->getPosition();
				}

				std::vector<double> mutant = xr4;
				std::vector<double> subtraction(this->D);
				subtract(xr0, xr1, subtraction);
				scale(subtraction, Fs[i]);
				add(mutant, subtraction, mutant);
				subtract(xr2, xr3, subtraction);
				scale(subtraction, Fs[i]);
				add(mutant, subtraction, mutant);
				T* m = new T(mutant);

				mutants.push_back(m);
			}
			
			return mutants;
		}
};
