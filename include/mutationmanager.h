#pragma once
#include <algorithm>
#include <iostream>
#include "particle.h"
#include "constrainthandler.h"
#include <functional>
#include <limits>
#include "rng.h"
#include "util.h"

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
};

class MutationManager {
	protected:
		int const D;
		ConstraintHandler* deCH;
		std::vector<Particle*> genomes;
		std::vector<double> Fs;
		Particle* best;
		Particle* pBest;
		virtual Particle* mutate(int const i)=0;
	public:
		MutationManager(int const D, ConstraintHandler* deCH);
		virtual ~MutationManager();
		std::vector<Particle*> mutate(std::vector<Particle*>const& genomes, std::vector<double>& Fs);
		static MutationManager* createMutationManager(MutationType const mutationType, int const D, ConstraintHandler* deCH);
};

class Rand1MutationManager : public MutationManager {
	public:
		Rand1MutationManager(int const D, ConstraintHandler* deCH);
		Particle* mutate(int const i);
};

class TTB1MutationManager : public MutationManager {
	public:
		TTB1MutationManager(int const D, ConstraintHandler* deCH);
		Particle* mutate(int const i);
};

class TTPB1MutationManager : public MutationManager {
	public:
		TTPB1MutationManager(int const D, ConstraintHandler* deCH);
		Particle* mutate(int const i);
};


class Best1MutationManager: public MutationManager {
	public:
		Best1MutationManager(int const D, ConstraintHandler* deCH);
		Particle* mutate(int const i);
};

class Best2MutationManager: public MutationManager {
	public:
		Best2MutationManager(int const D, ConstraintHandler* deCH);
		Particle* mutate(int const i);
};

class Rand2MutationManager: public MutationManager {
	public:
		Rand2MutationManager(int const D, ConstraintHandler* deCH);
		Particle* mutate(int const i);
};

class Rand2DirMutationManager : public MutationManager {
	public:
		Rand2DirMutationManager(int const D, ConstraintHandler* deCH);
		Particle* mutate(int const i);
};

class NSDEMutationManager : public MutationManager {
	public:
		NSDEMutationManager(int const D, ConstraintHandler* deCH);
		Particle* mutate(int const i);
};

class TrigonometricMutationManager : public MutationManager {
	private:
		double const gamma;
		Particle* trigonometricMutation(int const i);
		Particle* rand1Mutation(int const i);

	public:
		TrigonometricMutationManager(int const D, ConstraintHandler* deCH);
		Particle* mutate(int const i);
};

class TwoOpt1MutationManager : public MutationManager {
	public:
		TwoOpt1MutationManager(int const D, ConstraintHandler* deCH);
		Particle* mutate(int const i);
};

class TwoOpt2MutationManager : public MutationManager {
	public:
		TwoOpt2MutationManager(int const D, ConstraintHandler* deCH);
		Particle* mutate(int const i);
};
