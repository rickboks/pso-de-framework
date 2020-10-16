#pragma once
#include<vector>
#include<functional>
#include<map>
#include<random>

class Particle;

class TopologyManager {
	protected:
		std::vector<Particle*>const &particles;
	public:
		TopologyManager(std::vector<Particle*> const & particles);
		virtual ~TopologyManager();
		virtual void update(double progress);
};

extern std::map<std::string, std::function<TopologyManager* (std::vector<Particle*> const&)>> const topologies;

class LbestTopologyManager : public TopologyManager {
	public:
		LbestTopologyManager(std::vector<Particle*> const & particles);
};

class GbestTopologyManager : public TopologyManager {
	public:
		GbestTopologyManager(std::vector<Particle*> const & particles);
};

class RandomTopologyManager : public TopologyManager {
	private:
		int const connections; //Amount of neigbors every particle gets
	public:
		RandomTopologyManager(std::vector<Particle*> const & particles);
};

class VonNeumannTopologyManager : public TopologyManager {
	public:
		VonNeumannTopologyManager(std::vector<Particle*> const & particles);
};

class WheelTopologyManager : public TopologyManager {
	public:
		WheelTopologyManager(std::vector<Particle*> const & particles);
};

class IncreasingTopologyManager : public TopologyManager {
	private:
		int currentConnectivity;
		int maxConnectivity; // Maximum amount of neighbors for a particle (popSize-1)
		int minConnectivity;
	public:
		IncreasingTopologyManager(std::vector<Particle*> const & particles);
		void update(double progress);
};

class DecreasingTopologyManager : public TopologyManager {
	private:
		int currentConnectivity;
		int maxConnectivity;
	public:
		DecreasingTopologyManager(std::vector<Particle*> const & particles);
		void update(double progress);
};

class MultiSwarmTopologyManager : public TopologyManager {
	private:
		int const clusterSize;
		int count;
		void createClusters();
	public:
		MultiSwarmTopologyManager(std::vector<Particle*> const & particles);
		void update(double progress);
};

