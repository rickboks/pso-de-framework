#pragma once
#include <vector>
#include <map>
#include <random>

class ConstraintHandler;
struct ParticleUpdateSettings;
class Particle;

enum UpdateManagerType {
	INERTIA_WEIGHT,
	DECR_INERTIA_WEIGHT,
	CONSTRICTION_COEFFICIENT,
	FIPS,
	BARE_BONES,
	MAN_END
};

class ParticleUpdateManager {
	protected:
		std::vector<double> & x;
		std::vector<double> & v;
		std::vector<double> const& p;
		std::vector<double> const& g;
		int const D;
	public:
		ParticleUpdateManager(std::vector<double>& x, std::vector<double>& v,
			std::vector<double>const& p, std::vector<double>const& g);
		virtual ~ParticleUpdateManager();

		//This constructor is used for ParticleUpdateManagers that do not use a 
		//velocity vector.
		ParticleUpdateManager(std::vector<double>& x, std::vector<double> & p, std::vector<double>& g);

		static ParticleUpdateManager* createParticleUpdateManager(std::vector<double>& x, 
			std::vector<double>& v,	std::vector<double>const & p, std::vector<double>const& g, 
			ParticleUpdateSettings const settings, std::vector<Particle*>& neighborhood);

		virtual void updateVelocity(double const progress);
		virtual void updatePosition();
};

class InertiaWeightManager : public ParticleUpdateManager{
	private:
		double const phi1;
		double const phi2;
		double w;

	public:
		InertiaWeightManager(std::vector<double>& x, std::vector<double>& v,
			std::vector<double>const& p, std::vector<double>const& g, std::map<int, double> paramaters);
		void updateVelocity(double const progress);
};

class DecrInertiaWeightManager : public ParticleUpdateManager {
	private:
		double const phi1;
		double const phi2;
		double const wMin;
		double const wMax;
	public:
		DecrInertiaWeightManager(std::vector<double>& x, std::vector<double>& v,
			std::vector<double>const& p, std::vector<double>const& g, std::map<int, double> paramaters);
		void updateVelocity(double const progress);

};

class ConstrictionCoefficientManager : public ParticleUpdateManager {
	private:
		double const phi1;
		double const phi2;
		double const chi;
	public: 

		ConstrictionCoefficientManager(std::vector<double> & x, std::vector<double> & v,
			std::vector<double>const& p, std::vector<double>const& g, std::map<int, double> paramaters);

		void updateVelocity(double const progress);
};

class FIPSManager : public ParticleUpdateManager {
	private:
		double const phi;
		double const chi;
		std::vector<Particle*>& neighborhood;
	public:
		FIPSManager(std::vector<double> & x, std::vector<double> & v,
			std::vector<double>const& p, std::vector<double>const& g, 
			std::map<int, double> paramaters, std::vector<Particle*>& neighborhood);
		void updateVelocity(double const progress);
};



class BareBonesManager : public ParticleUpdateManager {
	private:
		
	public: 
		BareBonesManager(std::vector<double> & x, std::vector<double> & v,
			std::vector<double>const& p, std::vector<double>const& g, std::map<int, double> paramaters);
		void updatePosition();
		void updateVelocity(double const progress);
};
