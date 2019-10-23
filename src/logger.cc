#include "logger.h"
#include <string>

Logger::Logger(std::vector<Solution*>& solutions)
: solutions(solutions){
	//outfile->open(outputFile);
}

void Logger::logPositions(){
	for (unsigned int i = 0; i < solutions.size(); i++){
		*outfile << solutions[i]->positionString() << std::endl;
	}

	*outfile << std::endl;
}

void Logger::logInfo(Problem* problem){
	*outfile << std::to_string(problem->smallestValues[0]) << " " << std::to_string(problem->largestValues[0]) << std::endl;
	*outfile << std::to_string(problem->smallestValues[1]) << " " << std::to_string(problem->largestValues[1]) << std::endl;
	*outfile << coco_problem_get_name(problem->PROBLEM) << std::endl;
}