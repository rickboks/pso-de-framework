#include "util.h"
#include "rng.h"
#include <experimental/filesystem>

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


bool comparePtrs(Solution*a, Solution*b){
	return *a < *b;
}

Solution* getPBest(std::vector<Solution*> genomes, double const p){
	sortOnFitness(genomes);
	std::vector<Solution*> bestP = std::vector<Solution*>(genomes.begin(), genomes.begin() + (genomes.size() * p));
	if (bestP.empty())
		return genomes[0];
	else
		return bestP[rng.randInt(0, bestP.size()-1)];
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

void printVec(std::vector<double> v){
	for (double d : v)
		std::cout << d << " ";
	std::cout << std::endl;
}

