#pragma once
#include <vector>
#include <set>
#include <map>
#include <functional>
#include <stdexcept>

class DEAdaptationManager {
protected:
	std::set<double> SCr;
	std::set<double> SF;

	std::vector<double> previousFs;
	std::vector<double> previousCrs;
public:
	virtual ~DEAdaptationManager(){};
	virtual void successfulIndex(int const i)=0;
	virtual void nextF(std::vector<double>& Fs)=0;
	virtual void nextCr(std::vector<double>& Crs)=0;
	virtual void reset()=0;
	virtual void update()=0;
};

extern std::map<std::string, std::function<DEAdaptationManager*()>> const deAdaptations;

class JADEManager : public DEAdaptationManager{
private:
	double MuCr;
	double MuF;
	double const c;
	double lehmerMean() const;
public:
	JADEManager();
	void successfulIndex(int const i);
	void nextF(std::vector<double>& Fs);
	void nextCr(std::vector<double>& Crs);
	void reset();
	void update();
};

class NoAdaptationManager : public DEAdaptationManager {
private:
	double const F;
	double const Cr;
public:
	NoAdaptationManager();
	void successfulIndex(int const i);
	void nextF(std::vector<double>& Fs);
	void nextCr(std::vector<double>& Crs);
	void reset();
	void update();
};
