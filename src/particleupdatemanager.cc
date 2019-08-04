#include <vector>
#include <random>
#include <algorithm>
#include <functional>
#include <iostream>
#include <math.h>
#include "particle.h"
#include "particleupdatemanager.h"
#include "particleupdatesettings.h"

/*		Base 		*/
ParticleUpdateManager::ParticleUpdateManager(std::vector<double>& x, std::vector<double>& v,
	std::vector<double> & p, std::vector<double>& g, std::vector<double> const vMax)
	:x(x), v(v), p(p), g(g), D(x.size()), vMax(vMax), useVMax(true), 
	randDev(), 
	generator(randDev()){
}

ParticleUpdateManager::~ParticleUpdateManager(){

}

void ParticleUpdateManager::applyVMax(){
	for (int i = 0; i < D; i++){
		if (v[i] > vMax[i])
			v[i] = vMax[i];
		else if (v[i] < -vMax[i])
			v[i] = -vMax[i];
	}
}

void ParticleUpdateManager::updatePosition(){
		std::transform (x.begin(), x.end(),
				v.begin(), x.begin(),
		std::plus<double>());
}

/*		Factory 		*/
ParticleUpdateManager* ParticleUpdateManagerFactory::createParticleUpdateManager
	(std::vector<double>& x, std::vector<double>& v,std::vector<double> & p, 
	std::vector<double>& g,  ParticleUpdateSettings const settings, std::vector<Particle*>& neighborhood){

	switch(settings.managerType){
		case UpdateManagerType::INERTIA_WEIGHT:
			return new InertiaWeightManager(x,v,p,g,settings.parameters, settings.vMax);
		case UpdateManagerType::DECR_INERTIA_WEIGHT:
			return new DecrInertiaWeightManager(x,v,p,g,settings.parameters, settings.vMax);
		case UpdateManagerType::VMAX:
			return new VmaxManager(x,v,p,g,settings.parameters, settings.vMax);
		case UpdateManagerType::CONSTRICTION_COEFFICIENT:
			return new ConstrictionCoefficientManager(x,v,p,g,settings.parameters, settings.vMax);
		case UpdateManagerType::FIPS:
			return new FIPSManager(x,v,p,g,settings.parameters, neighborhood, settings.vMax);
		case UpdateManagerType::BARE_BONES:
			return new BareBonesManager(x,v,p,g,settings.parameters, settings.vMax);
		default: 
			throw std::invalid_argument("Error: Invalid particle update manager type");
		
	}
}

/*		VMax 		*/
VmaxManager::VmaxManager(std::vector<double> & x, std::vector<double> & v,
	std::vector<double> & p, std::vector<double> & g,  std::map<int, double> parameters, std::vector<double> const vMax)
	: ParticleUpdateManager(x,v,p,g, vMax),
	phi1 (parameters.find(Setting::S_VMAX_PHI1) != parameters.end() ? parameters[Setting::S_VMAX_PHI1] : VMAX_PHI1_DEFAULT),
	phi2 (parameters.find(Setting::S_VMAX_PHI2) != parameters.end() ? parameters[Setting::S_VMAX_PHI2] : VMAX_PHI2_DEFAULT),
	distr1(0,phi1), distr2(0,phi2){
}

void VmaxManager::update(double progress){
	double const r1 = distr1(generator);
	double const r2 = distr2(generator);

	std::vector<double> pMinx;
	pMinx.resize(D);
	std::transform( p.begin(), p.end(),
	                x.begin(), pMinx.begin(), 
	                std::minus<double>());

	std::vector<double> gMinx;
	gMinx.resize(D);
	std::transform( g.begin(), g.end(),
	                x.begin(), gMinx.begin(), 
	                std::minus<double>());


	std::transform(pMinx.begin(), pMinx.end(), pMinx.begin(),
           std::bind(std::multiplies<double>(), std::placeholders::_1, r1));

	std::transform(gMinx.begin(), gMinx.end(), gMinx.begin(),
       std::bind(std::multiplies<double>(), std::placeholders::_1, r2));

	std::transform (v.begin(), v.end(),
					pMinx.begin(), v.begin(),
					std::plus<double>());

	std::transform (v.begin(), v.end(),
					gMinx.begin(), v.begin(),
					std::plus<double>());
	if (useVMax)
		applyVMax();

	updatePosition();
}


/*		Inertia weight 		*/
InertiaWeightManager::InertiaWeightManager (std::vector<double>& x, std::vector<double>& v,
	std::vector<double> & p, std::vector<double>& g,  std::map<int, double> parameters, std::vector<double> const vMax)
	: ParticleUpdateManager(x,v,p,g, vMax),
	phi1 (parameters.find(Setting::S_INER_PHI1) != parameters.end() ? parameters[Setting::S_INER_PHI1] : INER_PHI1_DEFAULT),
	phi2 (parameters.find(Setting::S_INER_PHI2) != parameters.end() ? parameters[Setting::S_INER_PHI2] : INER_PHI2_DEFAULT),	
	w (parameters.find(Setting::S_INER_W) != parameters.end() ? parameters[Setting::S_INER_W] : INER_W_DEFAULT),
	distr1(0,phi1),
	distr2(0,phi2){
}

void InertiaWeightManager::update(double progress) {
	double const r1 = distr1(generator);
	double const r2 = distr2(generator);

	std::vector<double> pMinx(D);

	std::transform( p.begin(), p.end(),
	                x.begin(), pMinx.begin(), 
	                std::minus<double>());

	std::vector<double> gMinx(D);
	
	std::transform( g.begin(), g.end(),
	                x.begin(), gMinx.begin(), 
	                std::minus<double>());

	std::transform(pMinx.begin(), pMinx.end(), pMinx.begin(),
           std::bind(std::multiplies<double>(), std::placeholders::_1, r1));

	std::transform(gMinx.begin(), gMinx.end(), gMinx.begin(),
       std::bind(std::multiplies<double>(), std::placeholders::_1, r2));

	std::transform(v.begin(), v.end(), v.begin(),
               std::bind(std::multiplies<double>(), std::placeholders::_1, w));

	std::transform (v.begin(), v.end(),
					pMinx.begin(), v.begin(),
					std::plus<double>());

	std::transform (v.begin(), v.end(),
					gMinx.begin(), v.begin(),
					std::plus<double>());

	if (useVMax)
		applyVMax();

	updatePosition();
}

/*	Decreasing inertia weight manager */
DecrInertiaWeightManager::DecrInertiaWeightManager (std::vector<double>& x, std::vector<double>& v,
	std::vector<double> & p, std::vector<double>& g,  std::map<int, double> parameters, std::vector<double> const vMax)
	: ParticleUpdateManager(x,v,p,g, vMax),
	phi1 (parameters.find(Setting::S_DINER_PHI1) != parameters.end() ? parameters[Setting::S_DINER_PHI1] : DINER_PHI2_DEFAULT),
	phi2 (parameters.find(Setting::S_DINER_PHI2) != parameters.end() ? parameters[Setting::S_DINER_PHI2] : DINER_PHI2_DEFAULT),	
	w (parameters.find(Setting::S_DINER_W_START) != parameters.end() ? parameters[Setting::S_DINER_W_START] : DINER_W_START_DEFAULT),
	wMin (parameters.find(Setting::S_DINER_W_END) != parameters.end() ? parameters[Setting::S_DINER_W_END] : DINER_W_END_DEFAULT),
	wMax(w),
	distr1(0,phi1),
	distr2(0,phi2) {

}

void DecrInertiaWeightManager::update(double progress) {
	double const r1 = distr1(generator);
	double const r2 = distr2(generator);

	std::vector<double> pMinx(D);
	std::transform( p.begin(), p.end(),
	                x.begin(), pMinx.begin(), 
	                std::minus<double>());

	std::vector<double> gMinx(D);

	std::transform( g.begin(), g.end(),
	                x.begin(), gMinx.begin(), 
	                std::minus<double>());

	std::transform(pMinx.begin(), pMinx.end(), pMinx.begin(),
           std::bind(std::multiplies<double>(), std::placeholders::_1, r1));

	std::transform(gMinx.begin(), gMinx.end(), gMinx.begin(),
       std::bind(std::multiplies<double>(), std::placeholders::_1, r2));

	std::transform(v.begin(), v.end(), v.begin(),
               std::bind(std::multiplies<double>(), std::placeholders::_1, w));

	std::transform (v.begin(), v.end(),
					pMinx.begin(), v.begin(),
					std::plus<double>());

	std::transform (v.begin(), v.end(),
					gMinx.begin(), v.begin(),
					std::plus<double>());

	if (useVMax)
		applyVMax();

	w = wMax - progress * (wMax - wMin);
	std::cout << "w: " << w << std::endl;
	
	updatePosition();
}


/*		Constriction Coefficient 		*/
ConstrictionCoefficientManager::ConstrictionCoefficientManager(std::vector<double> & x, std::vector<double> & v,
	std::vector<double> & p, std::vector<double> & g,  std::map<int, double> parameters, std::vector<double> const vMax)
	: ParticleUpdateManager(x,v,p,g, vMax),
	phi1 (parameters.find(Setting::S_CC_PHI1) != parameters.end() ? parameters[Setting::S_CC_PHI1] : CC_PHI1_DEFAULT),
	phi2 (parameters.find(Setting::S_CC_PHI2) != parameters.end() ? parameters[Setting::S_CC_PHI2] : CC_PHI2_DEFAULT),
	chi (2 / ((phi1+phi2) -2 + sqrt( pow(phi1+phi2, 2) - 4 * (phi1+phi2)))),
	distr1(0,phi1),
	distr2(0,phi2){
}


void ConstrictionCoefficientManager::update(double progress){
	double const r1 = distr1(generator);
	double const r2 = distr2(generator);

	std::vector<double> pMinx(D);
	pMinx.resize(D);
	std::transform( p.begin(), p.end(),
	                x.begin(), pMinx.begin(), 
	                std::minus<double>());

	std::vector<double> gMinx(D);
	std::transform( g.begin(), g.end(),
	                x.begin(), gMinx.begin(), 
	                std::minus<double>());

	std::transform(pMinx.begin(), pMinx.end(), pMinx.begin(),
        std::bind(std::multiplies<double>(), std::placeholders::_1, r1));

	std::transform(gMinx.begin(), gMinx.end(), gMinx.begin(),
    	std::bind(std::multiplies<double>(), std::placeholders::_1, r2));

	std::transform (v.begin(), v.end(),
					pMinx.begin(), v.begin(),
					std::plus<double>());

	std::transform (v.begin(), v.end(),
					gMinx.begin(), v.begin(),
					std::plus<double>());

	std::transform(v.begin(), v.end(), v.begin(),
       std::bind(std::multiplies<double>(), std::placeholders::_1, chi));
	
	if (useVMax)
		applyVMax();

	updatePosition();
}

/*		Fully Informed 		*/
FIPSManager::FIPSManager(std::vector<double> & x, std::vector<double> & v,
	std::vector<double> & p, std::vector<double> & g,  std::map<int, double> parameters
	, std::vector<Particle*>& neighborhood, std::vector<double> const vMax)
	: ParticleUpdateManager(x,v,p,g, vMax),
	phi (parameters.find(Setting::S_FIPS_PHI) != parameters.end() ? parameters[Setting::S_FIPS_PHI] : FIPS_PHI_DEFAULT),
	chi (2 / ((phi) -2 + sqrt( pow(phi, 2) - 4 * (phi)))),
	neighborhood(neighborhood),
	distr(0,phi){
}


void FIPSManager::update(double progress){
	double const oneOverK = 1.0/neighborhood.size();

	std::vector<double> sum(D,0.0);
	std::vector<double> pMinx(D);

	std::vector< std::vector<double> > p_n;
	p_n.reserve(neighborhood.size());
	for (Particle* n : neighborhood)
		p_n.push_back(n->getP());

	for (int i = 0; i < (int) neighborhood.size(); i++){
		double const r = distr(generator);
		
		std::transform( p_n[i].begin(), p_n[i].end(),
		                x.begin(), pMinx.begin(), 
		                std::minus<double>());

		std::transform(pMinx.begin(), pMinx.end(), pMinx.begin(),
       		std::bind(std::multiplies<double>(), std::placeholders::_1, r));

		std::transform (sum.begin(), sum.end(),
						pMinx.begin(), sum.begin(),
						std::plus<double>());		
	}
	
	std::transform(sum.begin(), sum.end(), sum.begin(),
       		std::bind(std::multiplies<double>(), std::placeholders::_1, oneOverK));

	std::transform(v.begin(), v.end(),
					sum.begin(), v.begin(),
					std::plus<double>());

	std::transform(v.begin(), v.end(), v.begin(),
	std::bind(std::multiplies<double>(), std::placeholders::_1, chi));

	if (useVMax)
		applyVMax();

	updatePosition();
}

/* 		Bare Bones 		*/
BareBonesManager::BareBonesManager(std::vector<double> & x, std::vector<double> & v,
	std::vector<double> & p, std::vector<double> & g,  std::map<int, double> parameters,
	std::vector<double> const vMax) :
	ParticleUpdateManager(x,v,p,g,vMax) {

}

void BareBonesManager::update(double progress){
	for (int i = 0; i < D; i++){
		std::normal_distribution<double> distr((g[i] + p[i]) / 2.0, fabs(g[i] - p[i]));
		x[i] = distr(generator);
	}
}
