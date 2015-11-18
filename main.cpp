#include "netcart.h"
#include "util.h"
#include <iostream>


int main()
{
	Eigen::MatrixXd m;
	Netcart mynetcart("data/hospital.edges","data/hospital.nodefeat");
	mynetcart.initialization();
// give data to Rodan.py just for debug
	WriteMat(mynetcart.X, "/home/jackhaha363/netcart/netcartPython/X.txt");
	WriteMat(mynetcart.R, "/home/jackhaha363/netcart/netcartPython/R.txt");
	WriteMat(mynetcart.W, "/home/jackhaha363/netcart/netcartPython/W.txt");
	WriteMat(mynetcart.Constant, "/home/jackhaha363/netcart/netcartPython/Const.txt");	
//	cout << mynetcart.LogLikelihood() << endl;
	mynetcart.SetLearningRate(0.01,0.1,0.1);

	// test XGradient
	Eigen::MatrixXd Xgrad = Eigen::MatrixXd::Zero(mynetcart.X.rows(), mynetcart.X.cols());
	mynetcart.X_vGradient(Xgrad, 47);
	// test Wgradient
	Eigen::MatrixXd Wgrad = Eigen::MatrixXd::Zero(mynetcart.W.rows(), mynetcart.W.cols());
	Eigen::MatrixXd Constgrad = Eigen::MatrixXd::Zero(mynetcart.Constant.rows(), mynetcart.Constant.cols());
	mynetcart.W_bGradient(Wgrad, Constgrad, 0);
//	mynetcart.Optimize();
	return 0;	
}