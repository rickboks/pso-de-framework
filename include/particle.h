#pragma once
#include <vector>
#include <map>
#include "particleupdatesettings.h"
#include "IOHprofiler_experimenter.h"

class ParticleUpdateManager;

class Particle {
	private:		
		std::vector<double> x;
		std::vector<double> v;
		std::vector<double> p;		
		std::vector<double> g;

		double pbest;
		double gbest;
		bool evaluated;
		double fitness;

		std::vector<Particle*> neighborhood;
		ParticleUpdateManager* particleUpdateManager;
		ParticleUpdateSettings const * const settings;
		ConstraintHandler const * const psoCH;
	public:
		int const D;
		bool const isPSO;
		Particle(int const D, ParticleUpdateSettings const*const particleUpdateSettings);
		Particle(int const D);
		Particle(Particle const & other);
		Particle(std::vector<double> x); 
		~Particle();
		std::vector<double> getV() const;
		double getV(int const dim) const;
		void setV(std::vector<double> const v);
		void setV(int const dim, double val);
		void setX(std::vector<double> const x, double const fitness, bool const updateVelocity);
		void setX(std::vector<double> x);
		void addNeighbor(Particle* const neighbor);
		void removeNeighbor(Particle* const neighbor);
		void removeAllNeighbors();
		double getGbest() const;
		double getPbest() const;
		std::vector<double> getG() const;
		std::vector<double> getP() const;
		double getP(int const i) const;
		void updatePbest();
		void updateGbest();
		void updateVelocityAndPosition(double const progress);
		bool isNeighbor(Particle* particle) const;
		void replaceNeighbors(std::map<Particle*, Particle*> const mapping);
		int getAmountOfNeighbors();

		std::vector<double> getX() const;
		double evaluate (std::shared_ptr<IOHprofiler_problem<double> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger);
		double getFitness() const;
		void setFitness(double const d);
		std::string positionString() const;
		void randomize(std::vector<double> const lowerBounds, std::vector<double> const upperBounds);
		bool operator < (const Particle& s) const;
		void setX(int const dim, double const val);
		double getX(int const dim) const;

};
