#pragma once
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

struct result {
	vector<double> F;
	double smag;
};

/*Find the maximum value of a double vector
INPUT: vector<double>
OUTPUT: double 
*/
double max(vector<double>& v) {
	double output = v[0];
	for (int i = 1; i < v.size(); i++) {
		if (v[i] > output) {
			output = v[i];
		}
	}
	return output;
}

vector<double> absolute(const vector<double>& v) {
	vector<double> output(v.size());
	for (int i = 0; i < v.size(); i++) {
		if (v[i] < 0) {
			output[i] = -v[i];
		}
		else {
			output[i] = v[i];
		}
	}
	return output;
}

std::vector<double> addVectors(const std::vector<double>& a, const std::vector<double>& b) {
	if (a.size() != b.size()) {
		std::cout << "Error: vectors have different sizes" << std::endl;
		return std::vector<double>();
	}

	std::vector<double> result(a.size());
	for (int i = 0; i < a.size(); i++) {
		result[i] = a[i] + b[i];
	}
	return result;
}

std::vector<double> subtractVectors(const std::vector<double>& a, const std::vector<double>& b) {
	if (a.size() != b.size()) {
		std::cout << "Error: vectors have different sizes" << std::endl;
		return std::vector<double>();
	}

	std::vector<double> result(a.size());
	for (int i = 0; i < a.size(); i++) {
		result[i] = a[i] - b[i];
	}
	return result;
}