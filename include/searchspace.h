#pragma once
struct SearchSpace{
	SearchSpace(double const start, double const end)
		: start(start), end(end){};
	double const start;
	double const end;
};