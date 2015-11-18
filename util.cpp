#include "util.h"

Eigen::MatrixXd Sign(Eigen::MatrixXd &input)
{
	Eigen::MatrixXd result = input;
	for (int i = 0; i < input.rows(); ++i)
		for (int j = 0; j < input.cols(); ++j)
			result(i,j) = input(i,j) >= 0? 1 : -1;
	return result;
}

void normalization(Eigen::MatrixXd &input)
{
	double sum = 0;
	for (int i = 0; i < input.rows(); ++i)
	{
		for (int j = 0; j < input.cols(); ++j)
		{
			sum += pow(input(i,j), 2.);
		}
	}
	assert(sum >= 0);
	if (sum != 0)
		input /= sqrt(sum);
}

void normalization(Eigen::VectorXd &input)
{
	double sum = 0;
	for (int i = 0; i < input.rows(); ++i)
	{
		sum += pow(input(i), 2.);
	}
	assert(sum > 0);
	input /= sqrt(sum);
}

void WriteMat(const Eigen::MatrixXd &input, string filename)
{
	fstream output;
	output.open(filename, std::fstream::out);
	if (!output.is_open())
	{
		cout << "error open file " << filename << endl;
		return;
	} 
	for (int i = 0; i < input.rows(); ++i)
	{
		for (int j = 0; j < input.cols(); ++j)
		{
			output << input(i,j) << " ";
		}
		output << endl;
	}
	output.close();
}

void WriteMat(const Eigen::VectorXd &input, string filename)
{
	fstream output;
	output.open(filename, std::fstream::out);
	if (!output.is_open())
	{
		cout << "error open file " << filename << endl;
		return;
	} 
	for (int i = 0; i < input.rows(); ++i)
	{
			output << input(i) << " ";
	}
	output << endl;
	output.close();
}

double L1Norm(Eigen::MatrixXd mat)
{
	double result = 0;
	for (int i = 0; i < mat.rows(); ++i)
		for (int j = 0; j < mat.cols(); ++j)
			result += abs(mat(i,j));
	return result;
}