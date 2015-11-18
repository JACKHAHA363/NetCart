//This file contains function that would be used
#ifndef UTIL_H
#define UTIL_H

#include "netcart.h"

using std::string;

void normalization(Eigen::MatrixXd &input);

void normalization(Eigen::VectorXd &input);

void WriteMat(const Eigen::MatrixXd &input, string filename);

void WriteMat(const Eigen::VectorXd &input, string filename);

double L1Norm(Eigen::MatrixXd mat);

Eigen::MatrixXd Sign(Eigen::MatrixXd &input);	

#endif