#pragma once
#include <map>
#include <cstdlib>
#include <vector>
#include "particleupdatemanager.h"
#include "constrainthandler.h"

constexpr double INER_PHI1_DEFAULT = 1.49618;
constexpr double INER_PHI2_DEFAULT = 1.49618;
constexpr double DINER_PHI1_DEFAULT = 1.49618;
constexpr double DINER_PHI2_DEFAULT = 1.49618;
constexpr double INER_W_DEFAULT = 0.7298;
constexpr double DINER_W_DEFAULT = 0.7298;
constexpr double INER_ITERATIONS_DEFAULT = 0;
constexpr double DINER_W_START_DEFAULT = 0.9;
constexpr double DINER_W_END_DEFAULT = 0.4;
constexpr double VMAX_PHI1_DEFAULT = 1.49618;
constexpr double VMAX_PHI2_DEFAULT = 1.49618;
constexpr double CC_PHI1_DEFAULT = 2.05;
constexpr double CC_PHI2_DEFAULT = 2.05;
constexpr double FIPS_PHI_DEFAULT = 4.1;

enum Setting {
	S_INER_PHI1,
	S_INER_PHI2,
	S_DINER_PHI1,
	S_DINER_PHI2,
	S_VMAX_PHI1,
	S_VMAX_PHI2,
	S_CC_PHI1,
	S_CC_PHI2,
	S_INER_W, // Inertia weight
	S_DINER_W, // Decreasing nertia weight
	S_VMAX,
	S_FIPS_PHI,
	S_DINER_W_START, // Starting value of decreasing inertia weight
	S_DINER_W_END	// Ending value of decreasing inertia weight
};

struct ParticleUpdateSettings {
	ParticleUpdateSettings(std::string const managerType, std::map<int, double> const parameters, PSOConstraintHandler*const psoCH)
		:managerType(managerType), parameters(parameters), psoCH(psoCH){
	};

	ParticleUpdateSettings(){}

	std::string managerType;
	std::map<int, double> parameters;
	PSOConstraintHandler* psoCH;
};
