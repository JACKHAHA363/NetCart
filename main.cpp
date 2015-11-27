#include "netcart.h"
#include "util.h"
#include <iostream>


int main()
{
	Eigen::MatrixXd m;
	Netcart mynetcart("data/hospital.edges","data/hospital.nodefeat");
	mynetcart.initialization();
	mynetcart.SetRegularization();
// give data to Rodan.py just for debug
//	WriteMat(mynetcart.X, "/home/jackhaha363/netcart/netcartPython/X.txt");
//	WriteMat(mynetcart.R, "/home/jackhaha363/netcart/netcartPython/R.txt");
//	WriteMat(mynetcart.W, "/home/jackhaha363/netcart/netcartPython/W.txt");
//	WriteMat(mynetcart.Constant, "/home/jackhaha363/netcart/netcartPython/Const.txt");	
//	cout << mynetcart.LogLikelihoodGraph() << endl;

	// test XGradient
//	Eigen::MatrixXd Xgrad = Eigen::MatrixXd::Zero(mynetcart.X.rows(), 1);
//	mynetcart.X_vGradient(Xgrad, 47);
	
	// test Wgradient
//	Eigen::MatrixXd Wgrad = Eigen::MatrixXd::Zero(mynetcart.W.rows(), 1);
//	Eigen::MatrixXd Constgrad = Eigen::MatrixXd::Zero(mynetcart.Constant.rows(), 1);
//	mynetcart.W_bGradient(Wgrad, Constgrad, 0);

	// test likelihood
//	mynetcart.LogLikelihood();

	// test costfunction
//	cout << mynetcart.CostFunction() << endl;

	mynetcart.SetLearningRate(0.01,0.01,0.1);
	mynetcart.Optimize(200,200,200);
	return 0;	
}
