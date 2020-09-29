#pragma once
#include <vector>
#include <map>
#include "particleupdatesettings.h"
#include "IOHprofiler_experimenter.h"
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
		ParticleUpdateSettings const * const settings;
		PSOConstraintHandler* const psoCH;
	public:
		Particle(int const D, ParticleUpdateSettings const*const particleUpdateSettings);
		Particle(Particle const & other);
		~Particle();
		std::vector<double> getV() const;
		double getV(int const dim) const;
		void setV(std::vector<double> const v);
		void setV(int const dim, double val);
		void setXandUpdateV(std::vector<double> const x, double const fitness); 
		double getGbest() const;
		double getPbest() const;
		std::vector<double> getG() const;
		std::vector<double> getP() const;
		double getP(int const i) const;
		void updatePbest();
		void updateGbest();
		void updateVelocityAndPosition(double const progress);
		void addNeighbor(Particle* const neighbor);
		void removeNeighbor(Particle* const neighbor);
		void removeAllNeighbors();
		bool isNeighbor(Particle* particle) const;
		int getNumberOfNeighbors() const;
};
