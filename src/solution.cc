#include "iohsrc/Template/IOHprofiler_problem.hpp"
#include "iohsrc/Template/Loggers/IOHprofiler_csv_logger.h"
#include "solution.h"
#include "rng.h"
#include <iostream>
#include "particle.h"

Solution::Solution(int const D)
:x(D), D(D), evaluated(false) {

}

Solution::Solution(Particle* particle)
	:x(particle->getPosition()), D(particle->getPosition().size()),  evaluated(true), fitness(particle->getFitness()){

}

Solution::Solution(std::vector<double>x): 
	x(x), D(x.size()), evaluated(false){
}

Solution::Solution(std::vector<double> x, double fitness): 
	x(x), D(x.size()), evaluated(true), fitness(fitness){
}


double Solution::getFitness() const{
	return fitness;
}

double Solution::evaluate(std::shared_ptr<IOHprofiler_problem<double> > problem, std::shared_ptr<IOHprofiler_csv_logger> logger) {
	if (!evaluated){
		evaluated = true;		
		fitness = problem->evaluate(x);
		logger->write_line(problem->loggerInfo());
		return fitness;
	} else {
		return fitness;
	}
}

std::vector<double> Solution::getPosition() const {
	return x;
}

std::string Solution::positionString() const {
	std::string pos = "";
	for (int i = 0; i < D -1; i++){
		pos += std::to_string(x[i]);
		pos += " ";
	}
	pos += std::to_string(x[D-1]);

	return pos;
}

void Solution::setPosition(std::vector<double> position, double fitness){
	evaluated = true;
	this->fitness = fitness;
	this->x = position;
}

void Solution::randomize(std::vector<double> lowerBounds, std::vector<double> upperBounds){
	for (int i = 0; i < D; i++){
		x[i] = rng.randDouble(lowerBounds[i], upperBounds[i]);
	}
}

int Solution::getDimension(){
	return D;
}

bool Solution::operator < (const Solution& s) const {
	return fitness < s.getFitness();
}