#include "deadaptationmanager.h"
#include <algorithm>
#include <numeric>
#include <iostream>
#include "rng.h"

#define LC(X) [](){return new X();}
std::map<std::string, std::function<DEAdaptationManager*()>> const deAdaptations({
		{"jade", LC(JADEManager)},
		{"none", LC(NoAdaptationManager)},
});

//JADE
JADEManager::JADEManager()
 : DEAdaptationManager(), MuCr(0.5), MuF(0.6), c(0.1){}

void JADEManager::successfulIndex(int i){
	succesfulValues(previousFs[i], previousCrs[i]);
}

void JADEManager::succesfulValues(double F, double Cr){
	SF.insert(F); SCr.insert(Cr);
};

void JADEManager::update(){
	if (!SCr.empty()){
		MuCr = (1.0-c) * MuCr + c * (std::accumulate(SCr.begin(), SCr.end(), 0.0)/SCr.size());
		MuCr = std::max(MuCr, 0.01);
		MuCr = std::min(MuCr, 1.0);
	}

	if (!SF.empty()){
		MuF = (1.0-c) * MuF + c * lehmerMean();
		MuF = std::max(MuF, 0.01);
		MuF = std::min(MuF, 1.2);
	}
}

void JADEManager::nextF(std::vector<double>& Fs){
	unsigned int third = Fs.size()/3;
	std::vector<double> indices(Fs.size());
	std::iota(indices.begin(), indices.end(), 0);
	rng.shuffle(indices.begin(), indices.end());

	for (unsigned int i = 0; i < third; i++)
		Fs[indices[i]] = std::max(rng.randDouble(0.0,1.2),0.0);

	for (unsigned int i = third; i < Fs.size(); i++)
		Fs[indices[i]] = std::max(rng.normalDistribution(MuF, 0.1),0.0);	

	previousFs = Fs;
}


void JADEManager::nextCr(std::vector<double>& Crs){
	for (unsigned int i = 0; i < Crs.size(); i++){
		Crs[i] = std::max(rng.normalDistribution(MuCr, 0.1),0.0);
		
	}
	previousCrs = Crs;
}

double JADEManager::lehmerMean(){
	double sumOfSquares = std::inner_product(SF.begin(), SF.end(), SF.begin(), 0.0);
	double sum = std::accumulate(SF.begin(), SF.end(), 0.0);
	return sumOfSquares/sum;
}

void JADEManager::reset(){
	SCr.clear();
	SF.clear();
}

//NO ADAPTATION
NoAdaptationManager::NoAdaptationManager()
 : DEAdaptationManager(), F(0.5), Cr(0.9){}

void NoAdaptationManager::succesfulValues(double F, double Cr){
	//ignore
}

void NoAdaptationManager::update(){
	//ignore
}

void NoAdaptationManager::nextF(std::vector<double>& Fs){
	std::fill(Fs.begin(), Fs.end(), F);
}

void NoAdaptationManager::nextCr(std::vector<double>& Crs){
	std::fill(Crs.begin(), Crs.end(), Cr);
}

void NoAdaptationManager::reset(){
	//ignore
}

void NoAdaptationManager::successfulIndex(int i){
	//ignore
}
