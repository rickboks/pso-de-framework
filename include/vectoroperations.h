#pragma once
#include <vector>
#include <functional>
#include <algorithm>

void scale(std::vector<double>& vec, double x);
void add(std::vector<double>const& lhs, std::vector<double>const& rhs, std::vector<double>& store);
void subtract(std::vector<double>const& lhs, std::vector<double>const& rhs, std::vector<double>& store);
void randomMult(std::vector<double>& vec, double min, double max);
