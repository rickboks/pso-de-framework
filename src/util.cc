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

Particle* getBest(std::vector<Particle*>const& genomes){
	int best = 0;
	double bestF = std::numeric_limits<double>::max();

	for (unsigned int i = 0; i < genomes.size(); i++){
		double const score = genomes[i]->getFitness();
		if (score < bestF){
			bestF = score;
			best = i;
		}
	}
	return genomes[best];
}

Particle* getWorst(std::vector<Particle*>const& genomes){
	int worst = 0;
	double worstF = -std::numeric_limits<double>::max();

	for (unsigned int i = 0; i < genomes.size(); i++){
		double const score = genomes[i]->getFitness();
		if (score > worstF){
			worstF = score;
			worst = i;
		}
	}

	return genomes[worst];
}

Particle* pickRandom(std::vector<Particle*>& possibilities){
	int const r = rng.randInt(0,possibilities.size()-1);
	Particle* g = possibilities[r];
	possibilities.erase(possibilities.begin() + r);
	return g;
}

std::vector<Particle*> pickRandom(std::vector<Particle*>& possibilities, int n){
	std::vector<Particle*> particles;
	for (int i = 0; i < n; i++){
		particles.push_back(pickRandom(possibilities));
	}
	return particles;
}


bool comparePtrs(Particle*a, Particle*b){
	return *a < *b;
}

void sortOnFitness(std::vector<Particle*>& genomes){
	std::sort(genomes.begin(), genomes.end(), comparePtrs);
}

Particle* getPBest(std::vector<Particle*> genomes, double const p){
	sortOnFitness(genomes);
	std::vector<Particle*> bestP = std::vector<Particle*>(genomes.begin(), genomes.begin() + (genomes.size() * p));
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
