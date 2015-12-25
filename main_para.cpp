#include "netcart_para.h"
#include "util.h"
#include <iostream>


int main()
{
	Eigen::MatrixXd m;
	Netcart mynetcart("data/enron.edges","data/enron.nodefeat");
	mynetcart.initialization();
	mynetcart.SetRegularization();

	/*
	// test XGradient
	Eigen::MatrixXd Xgrad = Eigen::MatrixXd::Zero(mynetcart.X.rows(), 1);
	mynetcart.X_vGradient(Xgrad, 47);
	
	// test Wgradient
	Eigen::MatrixXd Wgrad = Eigen::MatrixXd::Zero(mynetcart.W.rows(), 1);
	Eigen::MatrixXd Constgrad = Eigen::MatrixXd::Zero(mynetcart.Constant.rows(), 1);
	mynetcart.W_bGradient(Wgrad, Constgrad, 0);

	// test likelihood
	mynetcart.LogLikelihood();
    */
	// test costfunction
    for (int i = 0; i < 2000; i++)
        mynetcart.CostFunction();
    cout << mynetcart.CostFunction() << endl;
	

	// These parameters are set to be consistent with what that is in rodan.py
	//mynetcart.SetLearningRate(0.01,0.01,0.1);
	//mynetcart.Optimize(200,200,200);
	return 0;	
}
