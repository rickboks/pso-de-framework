#include <vector>
#include <random>
#include <algorithm>
#include <functional>
#include <iostream>
#include "particle.h"
#include "particleupdatemanager.h"
#include "particleupdatesettings.h"
#include "util.h"
#include "rng.h"

/*		Base 		*/
ParticleUpdateManager::ParticleUpdateManager(std::vector<double>& x, std::vector<double>& v,
	std::vector<double>const& p, std::vector<double>const& g)
	:x(x), v(v), p(p), g(g), D(x.size()){
}

ParticleUpdateManager::~ParticleUpdateManager(){}

void ParticleUpdateManager::updatePosition(){
		std::transform (x.begin(), x.end(),
				v.begin(), x.begin(),
		std::plus<double>());
}

void ParticleUpdateManager::updateVelocity(double const progress){
		std::transform (x.begin(), x.end(),
				v.begin(), x.begin(),
		std::plus<double>());
}

#define LC(X) [](std::vector<double>& x, std::vector<double>& v,\
			std::vector<double>const& p, std::vector<double>const& g, \
			std::map<int, double> parameters, std::vector<Particle*>& neighborhood){return new X(x,v,p,g, parameters, neighborhood);}

std::map<std::string, std::function<ParticleUpdateManager* (std::vector<double>&, std::vector<double>&,
		std::vector<double>const&, std::vector<double>const&, std::map<int,double>, std::vector<Particle*>&)>> const updateManagers({
		{"I", LC(InertiaWeightManager)},
		{"D", LC(DecrInertiaWeightManager)},
		{"C", LC(ConstrictionCoefficientManager)},
		{"F", LC(FIPSManager)},
		{"B", LC(BareBonesManager)}
});

/*		Inertia weight 		*/
InertiaWeightManager::InertiaWeightManager (std::vector<double>& x, std::vector<double>& v,
	std::vector<double>const& p, std::vector<double>const& g,  std::map<int, double> parameters, std::vector<Particle*>& neighborhood)
	: ParticleUpdateManager(x,v,p,g),
	phi1 (parameters.find(Setting::S_INER_PHI1) != parameters.end() ? parameters[Setting::S_INER_PHI1] : INER_PHI1_DEFAULT),
	phi2 (parameters.find(Setting::S_INER_PHI2) != parameters.end() ? parameters[Setting::S_INER_PHI2] : INER_PHI2_DEFAULT),	
	w (parameters.find(Setting::S_INER_W) != parameters.end() ? parameters[Setting::S_INER_W] : INER_W_DEFAULT){}

void InertiaWeightManager::updateVelocity(double const progress) {
	std::vector<double> pMinx(D);
	subtract(p,x,pMinx);
	std::vector<double> gMinx(D);	
	subtract(g,x,gMinx);
	randomMult(pMinx, 0, phi1);
	randomMult(gMinx, 0, phi2);	
	scale(v, w);
	add(v,pMinx, v);
	add(v,gMinx, v);
}

/*	Decreasing inertia weight manager */
DecrInertiaWeightManager::DecrInertiaWeightManager (std::vector<double>& x, std::vector<double>& v,
	std::vector<double>const& p, std::vector<double>const& g,  std::map<int, double> parameters, std::vector<Particle*>& neighborhood)
	: ParticleUpdateManager(x,v,p,g),
	phi1 (parameters.find(Setting::S_DINER_PHI1) != parameters.end() ? parameters[Setting::S_DINER_PHI1] : DINER_PHI2_DEFAULT),
	phi2 (parameters.find(Setting::S_DINER_PHI2) != parameters.end() ? parameters[Setting::S_DINER_PHI2] : DINER_PHI2_DEFAULT),	
	wMin (parameters.find(Setting::S_DINER_W_END) != parameters.end() ? parameters[Setting::S_DINER_W_END] : DINER_W_END_DEFAULT),
	wMax(parameters.find(Setting::S_DINER_W_START) != parameters.end() ? parameters[Setting::S_DINER_W_START] : DINER_W_START_DEFAULT){}

void DecrInertiaWeightManager::updateVelocity(double const progress) {
	std::vector<double> pMinx(D);
	subtract(p,x,pMinx);
	std::vector<double> gMinx(D);
	subtract(g,x,gMinx);
	randomMult(pMinx, 0, phi1);
	randomMult(gMinx, 0, phi2);
	scale(v,wMax - progress * (wMax - wMin));
	add(v,pMinx, v);
	add(v,gMinx,v);
}

/*		Constriction Coefficient 		*/
ConstrictionCoefficientManager::ConstrictionCoefficientManager(std::vector<double> & x, std::vector<double> & v,
	std::vector<double>const& p, std::vector<double>const& g,  std::map<int, double> parameters, std::vector<Particle*>& neighborhood)
	: ParticleUpdateManager(x,v,p,g),
	phi1 (parameters.find(Setting::S_CC_PHI1) != parameters.end() ? parameters[Setting::S_CC_PHI1] : CC_PHI1_DEFAULT),
	phi2 (parameters.find(Setting::S_CC_PHI2) != parameters.end() ? parameters[Setting::S_CC_PHI2] : CC_PHI2_DEFAULT),
	chi (2.0 / ((phi1+phi2) - 2 + sqrt(pow(phi1+phi2, 2.0) - 4 * (phi1+phi2)))){}

void ConstrictionCoefficientManager::updateVelocity(double const progress){
	std::vector<double> pMinx(D);
	subtract(p,x,pMinx);
	std::vector<double> gMinx(D);
	subtract(g,x,gMinx);
	randomMult(pMinx, 0, phi1);
	randomMult(gMinx, 0, phi2);
	add(v,pMinx,v);
	add(v,gMinx,v);
	scale(v, chi);
}

/*		Fully Informed 		*/
FIPSManager::FIPSManager(std::vector<double> & x, std::vector<double> & v,
	std::vector<double>const& p, std::vector<double>const& g,  std::map<int, double> parameters
	, std::vector<Particle*>& neighborhood)
	: ParticleUpdateManager(x,v,p,g),
	phi (parameters.find(Setting::S_FIPS_PHI) != parameters.end() ? parameters[Setting::S_FIPS_PHI] : FIPS_PHI_DEFAULT),
	chi (2.0 / ((phi) -2 + sqrt( pow(phi, 2.0) - 4 * (phi)))),
	neighborhood(neighborhood){}


void FIPSManager::updateVelocity(double const progress){
	std::vector<double> sum(D,0.0);
	std::vector<double> pMinx(D);

	std::vector< std::vector<double> > p_n;
	p_n.reserve(neighborhood.size());
	for (Particle* n : neighborhood)
		p_n.push_back(n->getP());

	for (int i = 0; i < (int) neighborhood.size(); i++){
		subtract(p_n[i], x, pMinx);
		scale(pMinx, rng.randDouble(0,phi));
		add(sum, pMinx, sum);		
	}

	scale(sum, 1.0/neighborhood.size());
	add(v,sum,v);
	scale(v, chi);
}

/* 		Bare Bones 		*/
BareBonesManager::BareBonesManager(std::vector<double> & x, std::vector<double> & v,
	std::vector<double>const& p, std::vector<double>const& g,  std::map<int, double> parameters, std::vector<Particle*>& neighborhood) :
	ParticleUpdateManager(x,v,p,g) {}

void BareBonesManager::updatePosition(){
	for (int i = 0; i < D; i++){
		x[i] = rng.normalDistribution((g[i] + p[i]) / 2.0, std::abs(g[i] - p[i]));
	}
}

void BareBonesManager::updateVelocity(double const progress){ /* Do nothing*/ }
