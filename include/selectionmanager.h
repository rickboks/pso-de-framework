#pragma once
#include "particle.h"
#include "deadaptationmanager.h"

enum SelectionType {
	P2, P3, U2, U3, SEL_END
};

class SelectionManager {
protected:
	int const D;
	DEAdaptationManager* dam;	
public:
	SelectionManager(int D, DEAdaptationManager* dam);
	void selection(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2);
	void checkSuccessfulIndices();
	virtual void select(std::vector<Particle*>& particles, 
		std::vector<Particle*>const& p0, std::vector<Particle*>const& p2) =0;

	static SelectionManager* createSelectionManager(SelectionType type, int const D, DEAdaptationManager* dam);
};


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
