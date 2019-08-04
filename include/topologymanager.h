#pragma once
#include<vector>
#include<random>

enum Topology {
	LBEST,
	GBEST,
	RANDOM_GRAPH,
	VON_NEUMANN,
	WHEEL,
	INCREASING,
	DECREASING,
	MULTI_SWARM,
	TOP_END
};

class Particle;


class TopologyManager {
	protected:
		std::vector<Particle*> const & particles;
	public:
		TopologyManager(std::vector<Particle*> const & particles);
		virtual ~TopologyManager();
		virtual void initialize() = 0;
		virtual void update(double progress);
		virtual int getClosestValidPopulationSize(int const popSize);
};

class TopologyManagerFactory{
	public:
		static TopologyManager* createTopologyManager(int const type, std::vector<Particle*> const & particles);
};

class LbestTopologyManager : public TopologyManager {
	public:
		LbestTopologyManager(std::vector<Particle*> const & particles);
		void initialize();
		int getClosestValidPopulationSize(int const popSize);
};

class GbestTopologyManager : public TopologyManager {
	public:
		GbestTopologyManager(std::vector<Particle*> const & particles);
		void initialize();
};

class RandomTopologyManager : public TopologyManager {
	private:
		int const connections; //Amount of neigbors every particle gets
	public:
		RandomTopologyManager(std::vector<Particle*> const & particles);
		void initialize();
		int getClosestValidPopulationSize(int const popSize);
};

class VonNeumannTopologyManager : public TopologyManager {
	public:
		VonNeumannTopologyManager(std::vector<Particle*> const & particles);
		void initialize();
};

class WheelTopologyManager : public TopologyManager {
	public:
		WheelTopologyManager(std::vector<Particle*> const & particles);
		void initialize();
};

class IncreasingTopologyManager : public TopologyManager {
	private:
		int currentConnectivity;
		int maxConnectivity; // Maximum amount of neighbors for a particle (popSize-1)
		int minConnectivity;
		std::random_device randDev;
		std::mt19937 generator;
	public:
		IncreasingTopologyManager(std::vector<Particle*> const & particles);
		void initialize();
		void update(double progress);
};

class DecreasingTopologyManager : public TopologyManager {
	private:
		int currentConnectivity;
		int maxConnectivity;
		std::random_device randDev;
		std::mt19937 generator;
	public:
		DecreasingTopologyManager(std::vector<Particle*> const & particles);
		void initialize();
		void update(double progress);
};

class MultiSwarmTopologyManager : public TopologyManager {
	private:
		int const clusterSize;
		int count;
		std::random_device randDev;
		std::mt19937 generator;
		void createClusters();
	public:
		MultiSwarmTopologyManager(std::vector<Particle*> const & particles);
		void initialize();
		void update(double progress);		;
		int getClosestValidPopulationSize(int const popSize);
};

