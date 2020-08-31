#pragma once
#include <vector>
#include <set>
#include <stdexcept>

enum DEAdaptationType {
	JADE, NO, DEA_END
};

class DEAdaptationManager {
protected:
	std::set<double> SCr;
	std::set<double> SF;

	std::vector<double> previousFs;
	std::vector<double> previousCrs;
public:
	virtual ~DEAdaptationManager(){};
	virtual void successfulIndex(int i)=0;
	virtual void succesfulValues(double F, double Cr)=0;
	virtual void nextF(std::vector<double>& Fs)=0;
	virtual void nextCr(std::vector<double>& Crs)=0;
	virtual void reset()=0;
	virtual void update()=0;
	static DEAdaptationManager* createDEAdaptationManager(DEAdaptationType type);
};

class JADEManager : public DEAdaptationManager{
private:
	double MuCr;
	double MuF;
	double c;
	double lehmerMean();
public:
	JADEManager();
	void successfulIndex(int i);
	void succesfulValues(double F, double Cr);
	void nextF(std::vector<double>& Fs);
	void nextCr(std::vector<double>& Crs);
	void reset();
	void update();
};

class NoAdaptationManager : public DEAdaptationManager {
private:
	double F;
	double Cr;
public:
	NoAdaptationManager();
	void successfulIndex(int i);
	void succesfulValues(double F, double Cr);
	void nextF(std::vector<double>& Fs);
	void nextCr(std::vector<double>& Crs);
	void reset();
	void update();
};
