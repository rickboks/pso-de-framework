#include <fstream>
#include <iostream>
#include <vector>
#include "util.h" 
class Solution;

class Logger { // Used to log arbitrary stuff
	private:
		std::ostream out;
	public:
		Logger(std::streambuf *const out): out(out){};
		~Logger(){};

		template <typename T>
		void log(std::vector<T*> const pop) {
			for (T* s : pop){
				for (double x : s->getX()){
					out << x << " ";
				}
				out << "\n";
			}
		}
};
