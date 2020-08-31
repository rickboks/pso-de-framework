#pragma once
#include <vector>
#include <map>
#include <particleupdatesettings.h>
#include "solution.h"

class ParticleUpdateManager;

class Particle : public Solution {
	private:		
		std::vector<double> v;
		std::vector<double> p;		
		std::vector<double> g;
		double pbest;
		double gbest;
		std::vector<Particle*> neighborhood;
		ParticleUpdateManager* particleUpdateManager;
		ParticleUpdateSettings settings;
		std::vector<double> vMax;		
	public:
		Particle(int const D, ParticleUpdateSettings& particleUpdateSettings);
		Particle(const Particle& other);
		Particle(std::vector<double> x); 
		~Particle();
		std::vector<double> getVelocity() const;
		void setVelocity(std::vector<double> v);
		void addNeighbor(Particle* const neighbor);
		void removeNeighbor(Particle* const neighbor);
		void removeAllNeighbors();
		double getGbest() const;
		double getPbest() const;
		std::vector<double> getG() const;
		std::vector<double> getP() const;
		void updatePbest();
		void updateGbest();
		void updateVelocityAndPosition(double progress);
		bool isNeighbor(Particle* particle) const;
		void replaceNeighbors(std::map<Particle*, Particle*> mapping);
		int getAmountOfNeighbors();
};
