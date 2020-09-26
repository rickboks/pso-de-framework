#pragma once
#include "particle.h"
#include "deadaptationmanager.h"

class SelectionManager {
protected:
	int const D;
	DEAdaptationManager* const dam;	
public:
	SelectionManager(int const D, DEAdaptationManager* const dam);
	void selection(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Solution*>const& p2) const;
	void checkSuccessfulIndices();
	virtual void select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Solution*>const& p2) const =0;
};

extern std::map<std::string, std::function<SelectionManager* (int const, DEAdaptationManager *)>> const selections;

class Pairwise2SelectionManager : public SelectionManager {
private:
public:
	Pairwise2SelectionManager(int D, DEAdaptationManager* dam);
	void select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Solution*>const& p2) const;
};

class Pairwise3SelectionManager : public SelectionManager {
private:
public:
	Pairwise3SelectionManager(int D, DEAdaptationManager* dam);
	void select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Solution*>const& p2) const;
};

class Union2SelectionManager : public SelectionManager {
private:
public:
	Union2SelectionManager(int D, DEAdaptationManager* dam);
	void select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Solution*>const& p2) const;
};

class Union3SelectionManager : public SelectionManager {
private:
public:
	Union3SelectionManager(int D, DEAdaptationManager* dam);
	void select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Solution*>const& p2) const;
};
