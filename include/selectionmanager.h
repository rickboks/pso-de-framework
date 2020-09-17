#pragma once
#include "particle.h"
#include "deadaptationmanager.h"

class SelectionManager {
protected:
	int const D;
	DEAdaptationManager* const dam;	
public:
	std::string shorthand;
	SelectionManager(int const D, DEAdaptationManager* const dam);
	void selection(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2);
	void checkSuccessfulIndices();
	virtual void select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2) =0;
};

extern std::map<std::string, std::function<SelectionManager* (int const, DEAdaptationManager *)>> const selections;

class Pairwise2SelectionManager : public SelectionManager {
private:
public:
	Pairwise2SelectionManager(int D, DEAdaptationManager* dam);
	void select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2);
};

class Pairwise3SelectionManager : public SelectionManager {
private:
public:
	Pairwise3SelectionManager(int D, DEAdaptationManager* dam);
	void select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2);
};

class Union2SelectionManager : public SelectionManager {
private:
public:
	Union2SelectionManager(int D, DEAdaptationManager* dam);
	void select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2);
};

class Union3SelectionManager : public SelectionManager {
private:
public:
	Union3SelectionManager(int D, DEAdaptationManager* dam);
	void select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2);
};
