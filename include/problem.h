#pragma once
#include "coco.h"
#include <limits>
#include <vector>
#include <iostream>
typedef void (*evaluate_function_t)(const double *x, double *y);

struct Problem {
	Problem(evaluate_function_t evalFunc, coco_problem_t * PROBLEM)
		: evalFunc(evalFunc), PROBLEM(PROBLEM){
			int const dimension = coco_problem_get_dimension(PROBLEM);
			double const* smallest = coco_problem_get_smallest_values_of_interest(PROBLEM);
			double const* largest = coco_problem_get_largest_values_of_interest(PROBLEM);

			smallestValues = std::vector<double>(smallest, smallest + dimension);
			largestValues = std::vector<double>(largest, largest + dimension);
		}; 

	evaluate_function_t evalFunc;
	coco_problem_t * PROBLEM;
	std::vector<double> smallestValues;
	std::vector<double> largestValues;
};