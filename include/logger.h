#include <fstream>
#include <string>
#include <vector>
#include "solution.h"
#include "problem.h"

class Logger {
	private:
		std::ofstream* outfile;
		std::vector<Solution*>& solutions;
	public:
		Logger();
		Logger(std::vector<Solution*>& solutions);
		void logPositions();
		void logInfo(Problem* problem);
};