#pragma once
#include <algorithm>
#include <iostream>
#include "particle.h"
#include "constrainthandler.h"
#include <functional>
#include <limits>
#include "rng.h"
#include "util.h"

class MutationManager {
	protected:
		int const D;
		DEConstraintHandler* const deCH;
		std::vector<Solution*> genomes;
		std::vector<double> Fs;
		Solution const* best;
		virtual Solution* mutate(int const i) const=0;
	public:
		MutationManager(int const D, DEConstraintHandler * const deCH):D(D), deCH(deCH){};
		virtual ~MutationManager(){};
		std::vector<Solution*> mutate(std::vector<Solution*>const& genomes, std::vector<double>const& Fs);
};

extern std::map<std::string, std::function<MutationManager* (int const, DEConstraintHandler* const)>> const mutations;

class Rand1MutationManager : public MutationManager {
	public:
		Rand1MutationManager(int const D, DEConstraintHandler* const deCH): MutationManager(D, deCH){};
		Solution* mutate(int const i) const;
};

class TTB1MutationManager : public MutationManager {
	public:
		TTB1MutationManager(int const D, DEConstraintHandler* const deCH): MutationManager(D, deCH){};
		Solution* mutate(int const i) const;
};

class TTPB1MutationManager : public MutationManager {
	public:
		TTPB1MutationManager(int const D, DEConstraintHandler* const deCH): MutationManager(D, deCH){};
		Solution* mutate(int const i) const;
};

class Best1MutationManager: public MutationManager {
	public:
		Best1MutationManager(int const D, DEConstraintHandler* const deCH):MutationManager(D, deCH){};
		Solution* mutate(int const i) const;
};

class Best2MutationManager: public MutationManager {
	public:
		Best2MutationManager(int const D, DEConstraintHandler* const deCH):MutationManager(D, deCH){};
		Solution* mutate(int const i) const;
};

class Rand2MutationManager: public MutationManager {
	public:
		Rand2MutationManager(int const D, DEConstraintHandler* const deCH):MutationManager(D, deCH){};
		Solution* mutate(int const i) const;
};

class Rand2DirMutationManager : public MutationManager {
	public:
		Rand2DirMutationManager(int const D, DEConstraintHandler* const deCH):MutationManager(D, deCH){};
		Solution* mutate(int const i) const;
};

class NSDEMutationManager : public MutationManager {
	public:
		NSDEMutationManager(int const D, DEConstraintHandler* const deCH):MutationManager(D, deCH){};
		Solution* mutate(int const i) const;
};

class TrigonometricMutationManager : public MutationManager {
	private:
		double const gamma;
		Solution* trigonometricMutation(int const i) const;
		Solution* rand1Mutation(int const i) const;

	public:
		TrigonometricMutationManager(int const D, DEConstraintHandler* const deCH): MutationManager(D, deCH), gamma(0.05){};
		Solution* mutate(int const i) const;
};

class TwoOpt1MutationManager : public MutationManager {
	public:
		TwoOpt1MutationManager(int const D, DEConstraintHandler* const deCH): MutationManager(D, deCH) {};
		Solution* mutate(int const i) const;
};

class TwoOpt2MutationManager : public MutationManager {
	public:
		TwoOpt2MutationManager(int const D, DEConstraintHandler* const deCH): MutationManager(D, deCH){};
		Solution* mutate(int const i) const;
};

class ProximityMutationManager : public MutationManager {
	public:
		ProximityMutationManager(int const D, DEConstraintHandler* const deCH): MutationManager(D, deCH){};
		Solution* mutate(int const i) const;
};
