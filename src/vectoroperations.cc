#include "vectoroperations.h"
#include "rng.h"

void scale(std::vector<double> & vec, double x){
	std::transform(vec.begin(), vec.end(), vec.begin(),
           std::bind(std::multiplies<double>(), std::placeholders::_1, x));
}

void add(std::vector<double>const& lhs, std::vector<double>const& rhs, std::vector<double>& store){
	std::transform( lhs.begin(), lhs.end(),
                rhs.begin(), store.begin(), 
                std::plus<double>());
}

void subtract(std::vector<double>const& lhs, std::vector<double>const& rhs, std::vector<double>& store){
	std::transform( lhs.begin(), lhs.end(),
	                rhs.begin(), store.begin(), 
	                std::minus<double>());
}

void randomMult(std::vector<double> & vec, double min, double max){
	for (unsigned int i = 0; i < vec.size(); i++){
		vec[i] *= rng.randDouble(min, max);
	}
}
