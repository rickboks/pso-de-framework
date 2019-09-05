#pragma once
#include <vector>
#include <functional>
#include <algorithm>

void multiply(std::vector<double> & vec, double x);
void add(std::vector<double> lhs, std::vector<double> rhs, std::vector<double>& store);
void subtract(std::vector<double> lhs, std::vector<double> rhs, std::vector<double>& store);
void randomMult(std::vector<double> & vec, double min, double max);
