#pragma once
#include "genome.h"
#include <random>

enum MutationType {
	RAND_1,
	BEST_1,
	TTB_1,
	BEST_2,
	RAND_2,
	RAND_2_DIR,
	TOPOLOGY,
	NSDE,
	MUT_END
};

class MutationManager;
class MutationManagerFactory {
	public:
		static MutationManager* createMutationManager(MutationType const mutationType, 
			std::vector<Genome*>& genomes, double const F, int const D);
};

class MutationManager {
	protected:
		std::vector<Genome*> & genomes;
		double const F;
		int const D;
		Genome* getBest();
		Genome* pickRandom(std::vector<Genome*>& possibilities);

	public:
		MutationManager(std::vector<Genome*>& genomes, double const F, int const D);
		virtual ~MutationManager();
		virtual std::vector<Genome*> mutate() = 0;
};

class Rand1MutationManager : public MutationManager {
	private:

	public:
		Rand1MutationManager(std::vector<Genome*>& genomes, double const F, int const D);
		std::vector<Genome*> mutate();
};

class TTB1MutationManager : public MutationManager {
	private:

	public:
		TTB1MutationManager(std::vector<Genome*>& genomes, double const F, int const D);
		std::vector<Genome*> mutate();
};

class Best1MutationManager: public MutationManager {
	private:

	public:
		Best1MutationManager(std::vector<Genome*>& genomes, double const F, int const D);
		std::vector<Genome*> mutate();
};

class Best2MutationManager: public MutationManager {
	private:

	public:
		Best2MutationManager(std::vector<Genome*>& genomes, double const F, int const D);
		std::vector<Genome*> mutate();
};

class Rand2MutationManager: public MutationManager {
	private:

	public:
		Rand2MutationManager(std::vector<Genome*>& genomes, double const F, int const D);
		std::vector<Genome*> mutate();
};

class Rand2DirMutationManager : public MutationManager {
	private:

	public:
		Rand2DirMutationManager(std::vector<Genome*>& genomes, double const F, int const D);
		std::vector<Genome*> mutate();
};

class NSDEMutationManager : public MutationManager {
	private:

	public:
		NSDEMutationManager(std::vector<Genome*>& genomes, double const F, int const D);
		std::vector<Genome*> mutate();

};

class TopologyMutationManager : public MutationManager {
	private:
		int const radius;
		double const alpha;
		double const beta;
		std::vector<Genome*> getNeighbors(int const i) const;
	public:
		TopologyMutationManager(std::vector<Genome*>& genomes, double const F, int const D);
		std::vector<Genome*> mutate();
};