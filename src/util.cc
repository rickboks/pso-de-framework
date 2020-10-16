#include "util.h"
#include "rng.h"
#include <experimental/filesystem>

void scale(std::vector<double> & vec, double const x){
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

void randomMult(std::vector<double>& vec, double const min, double const max){
	for (unsigned int i = 0; i < vec.size(); i++){
		vec[i] *= rng.randDouble(min, max);
	}
}


bool comparePtrs(Solution const* const a, Solution const *const b){
	return *a < *b;
}

std::string generateConfig(std::string const templateFile, std::string const name){
	std::string const folder = "configurations";
	std::string const cfgFile = folder + "/" + name + ".ini";

	if (!std::experimental::filesystem::exists(folder)){
		std::cerr << "Creating directory \"" << folder << "\"."<< std::endl;
		std::experimental::filesystem::create_directory(folder);
	}

	std::ifstream src(templateFile, std::ios::binary);
    std::ofstream dst(cfgFile, std::ios::binary);

	dst << src.rdbuf() 
		<< "result_folder = " + name << std::endl 
		<< "algorithm_name = " + name << std::endl;

	return cfgFile;
}

void printVec(std::vector<double> const v){
	for (double d : v)
		std::cout << d << " ";
	std::cout << std::endl;
}

double distance(Solution const*const s1, Solution const*const s2) {
	double d=0;
	for (int i = 0; i < s1->D; i++){
		d += std::pow(s1->getX(i) - s2->getX(i), 2.); 
	}
	return std::sqrt(d);
}

std::string checkFilename(std::string const fn){
	std::string newFn = fn;
	int i = 1;
	while (std::experimental::filesystem::exists(newFn)){
		newFn = fn + "." + std::to_string(i);
		i++;
	}
	return newFn;
}
