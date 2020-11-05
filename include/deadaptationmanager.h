#pragma once
#include "solution.h"
#include <vector>
#include <set>
#include <map>
#include <functional>
#include <stdexcept>

class DEAdaptationManager {
protected:
	int const popSize;
public:
	DEAdaptationManager(int const popSize); 
	virtual ~DEAdaptationManager(){};
	virtual void nextF(std::vector<double>& Fs)=0;
	virtual void nextCr(std::vector<double>& Crs)=0;
	virtual void update(std::vector<double>const& targets, std::vector<double>const& trials)=0;
};

extern std::map<std::string, std::function<DEAdaptationManager*(int const)>> const deAdaptations;

class JADEManager : public DEAdaptationManager{
private:
	std::vector<double> previousFs;
	std::vector<double> previousCrs;

	double MuCr;
	double MuF;
	double const c;
	double lehmerMean(std::vector<double>const& SF) const;
public:
	JADEManager(int const popSize);
	void nextF(std::vector<double>& Fs);
	void nextCr(std::vector<double>& Crs);
	void update(std::vector<double>const& orig, std::vector<double>const& trials);
};

class SHADEManager : public DEAdaptationManager {
	private:
		std::vector<double> previousFs;
		std::vector<double> previousCrs;

		int const H;
		std::vector<double> MCr;
		std::vector<double> MF;
		std::vector<int> r;

		int k;
		double weightedLehmerMean(std::vector<double>const& x, std::vector<double>const& w) const;
		double weightedMean(std::vector<double>const& x, std::vector<double>const& w) const;
		std::vector<double> w(std::vector<double>const& delta) const;
	public:
		SHADEManager(int const popSize);
		void nextF(std::vector<double>& Fs);
		void nextCr(std::vector<double>& Crs);
		void update(std::vector<double>const& orig, std::vector<double>const& trials);
};

class NoAdaptationManager : public DEAdaptationManager {
private:
	double const F;
	double const Cr;
public:
	NoAdaptationManager(int const popSize);
	void nextF(std::vector<double>& Fs);
	void nextCr(std::vector<double>& Crs);
	void update(std::vector<double>const& orig, std::vector<double>const& trials);
};
