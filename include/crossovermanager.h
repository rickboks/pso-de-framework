#pragma once
#include <random>
#include <stdexcept>
#include <vector>
#include<algorithm>
#include "solution.h"
#include "rng.h"


enum CrossoverType {
	BINOMIAL,
	EXPONENTIAL,
	CROSS_END
};

template <class T>
class CrossoverManager {
	protected:
		int const D;
	public:
		CrossoverManager(int const D): D(D){}
		virtual ~CrossoverManager(){}

		virtual std::vector<T*> crossover(std::vector<T*>& targets, 
			std::vector<T*>& donors, std::vector<double>& Crs) = 0;

		virtual std::vector<double> singleCrossover(std::vector<double> target, 
			std::vector<double> donor, double const Cr) =0;
};

template <class T>
class BinomialCrossoverManager : public CrossoverManager<T> {
	private:

	public:
		BinomialCrossoverManager(int const D) : CrossoverManager<T>(D){}
		std::vector<T*> crossover(std::vector<T*>& genomes, std::vector<T*>& mutants, std::vector<double>& Crs){
			std::vector<T*> trials;
			trials.reserve(genomes.size());
			std::vector<double> x(this->D);

			for (unsigned int i = 0; i < genomes.size(); i++){
				std::vector<double> mutantX = mutants[i]->getPosition();
				std::vector<double> parentX = genomes[i]->getPosition();
				trials.push_back(new T(singleCrossover(genomes[i]->getPosition(), mutants[i]->getPosition(), Crs[i])));
			}

			return trials;
		}

		std::vector<double> singleCrossover(std::vector<double> target, std::vector<double> donor, double const Cr){
			std::vector<double> x(this->D);
			int const jrand = rng.randInt(0,this->D-1);
			for (int j = 0; j < this->D; j++){
				if (j == jrand || rng.randDouble(0,1) < Cr){
					x[j] = donor[j];
				} else {
					x[j] = target[j]; 
				}
			}
			return x;
		}
};

template <class T>
class ExponentialCrossoverManager : public CrossoverManager<T> {
	private:

	public:
		ExponentialCrossoverManager(int const D): CrossoverManager<T>(D){}
		std::vector<T*> crossover(std::vector<T*>& genomes, std::vector<T*>& mutants, std::vector<double>& Crs){
			std::vector<T*> trials;
			for (unsigned int i = 0; i < genomes.size(); i++){
				trials.push_back(new T(singleCrossover(genomes[i]->getPosition(), mutants[i]->getPosition(), Crs[i])));
			}

			return trials;
		}
		std::vector<double> singleCrossover(std::vector<double> target, std::vector<double> donor, double const Cr){
			std::vector<double> x(this->D);
			std::vector<int> mutantIndices;
			int const n = rng.randInt(0,this->D-1);

			int L = 0;
			do {
				L++;
			} while (rng.randDouble(0,1) < Cr && L < this->D);

			int steps = 0;
			for (int i = n; steps < L; i++){
				mutantIndices.push_back(i%this->D);
				steps++;
			}

			for (int i = 0; i < this->D; i++) {
				if (std::find(mutantIndices.begin(), mutantIndices.end(), i) != mutantIndices.end()){
					x[i] = donor[i];
				} else {
					x[i] = target[i];
				}
			}

			return x;
		}
};

class CrossoverManagerFactory {
	public:
		template <class T>
		static CrossoverManager<T>* createCrossoverManager(CrossoverType const crossoverType, int const D){
			switch(crossoverType){
				case BINOMIAL:
					return new BinomialCrossoverManager<T>(D);
				case EXPONENTIAL:
					return new ExponentialCrossoverManager<T>(D);
				default:
					throw std::invalid_argument("Error: Invalid DE crossover type");
			}
		}
};