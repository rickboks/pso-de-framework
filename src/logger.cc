#include "logger.h"
#include <iomanip>
#include <numeric>

void Logger::log(int const function, int const D, std::vector<double> const percCorrected, 
		std::vector<double> const bestX, double const bestF, int const numEvals){
	out.precision(3);

	out << function << " " << D << " ";
	for (double d : percCorrected)
		out << d << " ";
	for (double d : bestX)
		out << d << " ";
	out << bestF << " " << numEvals << "\n";
}

void Logger::start(int const function, int const D){
	out << function << " " << D << ":";
}

void Logger::log(std::vector<double> F, std::vector<double> Cr){
	out.precision(3);
	double const avgF = std::accumulate(F.begin(), F.end(), 0.) / double(F.size());
	double const avgCr = std::accumulate(Cr.begin(), Cr.end(), 0.) / double(Cr.size());
	out << avgF << " " << avgCr << ","; 
}

void Logger::newLine(){
	out << "\n";
}

Logger::~Logger(){
	out.close();
};
