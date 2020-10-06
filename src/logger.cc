#include "logger.h"

void Logger::log(int const function, int const D, std::vector<double> const percCorrected, 
		std::vector<double> const bestX, double const bestF, int const numEvals){

	out << function << " " << D << " ";
	for (double d : percCorrected)
		out << d << " ";
	for (double d : bestX)
		out << d << " ";
	out << bestF << " " << numEvals << "\n";
}
