#pragma once
#include <vector>
#include <iostream>
void printVector(std::vector<double> x){
	std::cout << "{";
	for (unsigned int i = 0; i < x.size() -1; i++){
		std::cout << x[i] << ", ";
	}

	std::cout << x[x.size()-1] << "}" << std::endl;

}