#include "netcart.h"
#include "util.h"
#include <iostream>


int main()
{
	Eigen::MatrixXd m;
	Netcart mynetcart("data/hospital.edges","data/hospital.nodefeat");
	mynetcart.initialization();
	mynetcart.SetRegularization();
	mynetcart.SetLearningRate(0.01,0.01,0.1);
	mynetcart.Optimize(200,200,200);
	return 0;	
}
