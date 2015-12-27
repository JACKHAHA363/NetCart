//This file contains function that would be used
#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include "Eigen/Dense"
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <assert.h>
#include <cfloat>
#include <set>

using std::string;
using std::cout;
using std::endl;
using std::fstream;
using std::abs;

void normalization(Eigen::MatrixXd &input);

void normalization(Eigen::VectorXd &input);

void WriteMat(const Eigen::MatrixXd &input, string filename);

void WriteMat(const Eigen::VectorXd &input, string filename);

double L1Norm(Eigen::MatrixXd mat);

Eigen::MatrixXd Sign(Eigen::MatrixXd &input);	

#endif
