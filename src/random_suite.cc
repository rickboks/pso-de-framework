#include "random_suite.h"
#include "util.h"

Random_suite::Random_suite(){
	std::vector<int> problem_id = {1,2,3};
	std::vector<int> instance_id = {1};
	std::vector<int> dimension = {30};

	IOHprofiler_set_suite_problem_id(problem_id);
	IOHprofiler_set_suite_instance_id(instance_id);
	IOHprofiler_set_suite_dimension(dimension);
	IOHprofiler_set_suite_name("random");
	registerProblem();
	loadProblem();
}

Random_suite::Random_suite(std::vector<int> problem_id, std::vector<int> instance_id, std::vector<int> dimension){
	IOHprofiler_set_suite_problem_id(problem_id);
	IOHprofiler_set_suite_instance_id(instance_id);
	IOHprofiler_set_suite_dimension(dimension);
	IOHprofiler_set_suite_name("random");
	registerProblem();
	loadProblem();
}

void Random_suite::registerProblem(){
	registerInFactory<IOHprofiler_problem<double>, f0> regf0("f0");
	registerInFactory<IOHprofiler_problem<double>, g0> regg0("g0");
	registerInFactory<IOHprofiler_problem<double>, h0> regh0("h0");

	mapIDTOName(1, "f0");
	mapIDTOName(2, "g0");
	mapIDTOName(3, "h0");
}

Random_suite* Random_suite::createInstance(){
	return new Random_suite();
}

Random_suite* Random_suite::createInstance(std::vector<int> problem_id, std::vector<int> instance_id, 
		std::vector<int> dimension){
	return new Random_suite(problem_id, instance_id, dimension);
}

f0::f0(int instance_id, int dimension) : rng(dev()), dist(0.,1.){
	IOHprofiler_set_instance_id(instance_id);
	IOHprofiler_set_problem_name("f0");
	IOHprofiler_set_number_of_objectives(1);
	IOHprofiler_set_number_of_variables(dimension);
	IOHprofiler_set_lowerbound(std::vector<double>(dimension,0.));
	IOHprofiler_set_upperbound(std::vector<double>(dimension,1.));
}

double f0::internal_evaluate(std::vector<double> const& x) {
	return dist(rng);
};

f0* f0::createInstance(int instance_id = 1, int dimension = 30){
	return new f0(instance_id, dimension);
}

g0::g0(int instance_id, int dimension) : rng(dev()), dist(0.,1.){
	IOHprofiler_set_instance_id(instance_id);
	IOHprofiler_set_problem_name("g0");
	IOHprofiler_set_number_of_objectives(1);
	IOHprofiler_set_number_of_variables(dimension);
	IOHprofiler_set_lowerbound(std::vector<double>(dimension, 0.));
	IOHprofiler_set_upperbound(std::vector<double>(dimension, 100.));
}

double g0::internal_evaluate(std::vector<double> const& x) {
	return dist(rng);
};

g0* g0::createInstance(int instance_id = 1, int dimension = 30){
	return new g0(instance_id, dimension);
}

h0::h0(int instance_id, int dimension) : rng(dev()), dist(0.,1.){
	IOHprofiler_set_instance_id(instance_id);
	IOHprofiler_set_problem_name("h0");
	IOHprofiler_set_number_of_objectives(1);
	IOHprofiler_set_number_of_variables(dimension);
	IOHprofiler_set_lowerbound(std::vector<double>(dimension, -0.6));
	IOHprofiler_set_upperbound(std::vector<double>(dimension, 0.4));
}

double h0::internal_evaluate(std::vector<double> const& x) {
	return dist(rng);
};

h0* h0::createInstance(int instance_id = 1, int dimension = 30){
	return new h0(instance_id, dimension);
}
