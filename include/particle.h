#pragma once
#include <vector>
#include <map>
#include <particleupdatesettings.h>
typedef void (*evaluate_function_t)(const double *x, double *y);

class ParticleUpdateManager;

class Particle {
	private:
		int const D;
		std::vector<double> x;
		std::vector<double> v;
		std::vector<double> p;
		double fitness;
		double pbest;
		std::vector<double> g;
		double gbest;
		std::vector<Particle*> neighborhood;
		ParticleUpdateManager* particleUpdateManager;
		ParticleUpdateSettings const settings;
		std::vector<double> const vMax;
		bool evaluated;
	public:

		Particle(int const D, ParticleUpdateSettings& particleUpdateSettings);
		Particle(const Particle& other);
		~Particle();
		void randomize(std::vector<double> lowerBounds, std::vector<double> upperBounds);
		std::vector<double> getPosition() const;
		void setPosition(std::vector<double> position, double fitness);
		void addNeighbor(Particle* const neighbor);
		void removeNeighbor(Particle* const neighbor);
		void removeAllNeighbors();
		double getGbest() const;
		double getPbest() const;
		double getFitness() const;
		std::vector<double> getG() const;
		std::vector<double> getP() const;
		void updatePbest();
		void updateGbest();
		void updateVelocityAndPosition(double progress);
		bool isNeighbor(Particle* particle) const;
		double getResultingVelocity() const;
		double evaluate(evaluate_function_t evalFunc);
		void replaceNeighbors(std::map<Particle*, Particle*> mapping);

};