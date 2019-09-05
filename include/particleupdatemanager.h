#pragma once
#include <vector>
#include <map>
#include <random>

class ParticleUpdateSettings;
class Particle;
enum UpdateManagerType {
	INERTIA_WEIGHT,
	DECR_INERTIA_WEIGHT,
	VMAX,
	CONSTRICTION_COEFFICIENT,
	FIPS,
	BARE_BONES,
	MAN_END
};

class ParticleUpdateManager {
	protected:
		std::vector<double> & x;
		std::vector<double> & v;
		std::vector<double> const & p;
		std::vector<double> const & g;
		int const D;
		std::vector<double> const vMax;
		bool const useVMax;
		void applyVMax();
		void updatePosition();
	public:
		ParticleUpdateManager(std::vector<double>& x, std::vector<double>& v,
			std::vector<double> & p, std::vector<double>& g, std::vector<double> const vMax);
		virtual ~ParticleUpdateManager();

		//This constructor is used for ParticleUpdateManagers that do not use a 
		//velocity vector.
		ParticleUpdateManager(std::vector<double>& x, std::vector<double> & p, std::vector<double>& g, std::vector<double> const vMax);
		virtual void update(double progress) = 0;
};

class ParticleUpdateManagerFactory {
	public:
		static ParticleUpdateManager* createParticleUpdateManager(std::vector<double>& x, 
			std::vector<double>& v,	std::vector<double> & p, std::vector<double>& g, 
			ParticleUpdateSettings const settings, std::vector<Particle*>& neighborhood);

};

class InertiaWeightManager : public ParticleUpdateManager{
	private:
		double const phi1;
		double const phi2;
		double w;

	public:
		InertiaWeightManager(std::vector<double>& x, std::vector<double>& v,
			std::vector<double> & p, std::vector<double>& g, std::map<int, double> paramaters, std::vector<double> const vMax);
		void update (double progress);
};

class DecrInertiaWeightManager : public ParticleUpdateManager {
	private:
		double const phi1;
		double const phi2;
		double w;
		double const wMin;
		double const wMax;

	public:
		DecrInertiaWeightManager(std::vector<double>& x, std::vector<double>& v,
			std::vector<double> & p, std::vector<double>& g, std::map<int, double> paramaters, std::vector<double> const vMax);
		void update (double progress);

};

class VmaxManager : public ParticleUpdateManager {
	private:		
		double const phi1;
		double const phi2;
	public:

		VmaxManager(std::vector<double> & x, std::vector<double> & v,
			std::vector<double> & p, std::vector<double> & g, 
			std::map<int, double> paramaters, std::vector<double> const vMax);

		void update(double progress);
};

class ConstrictionCoefficientManager : public ParticleUpdateManager {
	private:
		double const phi1;
		double const phi2;
		double const chi;
	public: 

		ConstrictionCoefficientManager(std::vector<double> & x, std::vector<double> & v,
			std::vector<double> & p, std::vector<double> & g, 
			std::map<int, double> paramaters, std::vector<double> const vMax);

		void update(double progress);
};

class FIPSManager : public ParticleUpdateManager {
	private:
		double const phi;
		double const chi;
		std::vector<Particle*>& neighborhood;
	public:
		FIPSManager(std::vector<double> & x, std::vector<double> & v,
			std::vector<double> & p, std::vector<double> & g, 
			std::map<int, double> paramaters, std::vector<Particle*>& neighborhood, std::vector<double> const vMax);
		void update(double progress);
};



class BareBonesManager : public ParticleUpdateManager {
	private:
		
	public: 
		BareBonesManager(std::vector<double> & x, std::vector<double> & v,
			std::vector<double> & p, std::vector<double> & g, 
			std::map<int, double> paramaters, std::vector<double> const vMax);
		void update(double progress);
};