#include <IOHprofiler_suite.h>
#include <random>

class Random_suite : public IOHprofiler_suite<double> {
	public:
		static Random_suite* createInstance();
		static Random_suite* createInstance(std::vector<int> problem_id, std::vector<int> instance_id, 
				std::vector<int> dimension);
		Random_suite(std::vector<int> problem_id, std::vector<int> instance_id, std::vector<int> dimension);
		Random_suite();
		void registerProblem();
};

class f0 : public IOHprofiler_problem<double> {
	private:
		std::random_device dev;
		std::mt19937 rng;
		std::uniform_real_distribution<double> dist;
	public:
		static f0 * createInstance(int instance_id, int dimension);
		f0(int instance_id, int dimension);
		double internal_evaluate(std::vector<double> const& x);
		~f0(){};
};

class g0 : public IOHprofiler_problem<double> {
	private:
		std::random_device dev;
		std::mt19937 rng;
		std::uniform_real_distribution<double> dist;
	public:
		static g0 * createInstance(int instance_id, int dimension);
		g0(int instance_id, int dimension);
		double internal_evaluate(std::vector<double> const& x);
		~g0(){};
};

class h0 : public IOHprofiler_problem<double> {
	private:
		std::random_device dev;
		std::mt19937 rng;
		std::uniform_real_distribution<double> dist;
	public:
		static h0 * createInstance(int instance_id, int dimension);
		h0(int instance_id, int dimension);
		double internal_evaluate(std::vector<double> const& x);
		~h0(){};
};
