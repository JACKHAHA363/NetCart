#ifndef NETCART_H
#define NETCART_H

// TODO: Problem is in optimizing X

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


const int MAX_ITERATE = 100;
const double ZERO = 1e-14;

using namespace std;
class Netcart
{
public:
	Netcart();

	Netcart(string graph_file, string attri_file, int _k_target = 4, 
			double _alpha = 0.5, bool _digraph = true);

	//set parameters
	void SetRegularization(double alpha_x = 0, double alpha_r = 0.9375, double alpha_w = 0);
	void SetAttributedGraph(string graph_file, string attri_file);
	void SetNumRole(int _k_target);
	void SetAttriWeight(double _alpha);
	void SetDigraph(bool _digraph);
	void SetLearningRate(double beta_x, double beta_r, double beta_w);

	void SetX(Eigen::MatrixXd _X){X = _X;}
	void SetR(Eigen::MatrixXd _R){R = _R;}
	void SetW(Eigen::MatrixXd _W){W = _W;}

	//getter function for debugging
	void PrintXRW();
	void PrintX();
	void PrintR();
	void PrintW();
	void SaveXRW(string path = "");//default save to current path
	void SaveX(string path = "");
	void SaveR(string path = "");
	void SaveW(string path = "");
	void SaveConst(string path = "");

	//optimizing
	void initialization();
	double LogLikelihood();
	double LogLikelihoodGraph();
	double LogLikelihoodAttri();
	double CostFunction();
	void Optimize(int max_iterate_x, int max_iterate_r, int max_iterate_w);
	void OptimizeR(int max_iterate_r);
	void OptimizeX(int max_iterate_x);
	void OptimizeW(int max_iterate_w);

//private:
	double predictEdge(int i, int j); // calculate the probability of edge from node i to node j
	double predictAttri(int v, int i); // calculate the probability of node v has attribute i
	bool Converge(double oldCost, double newCost); //test the converge condition
	void Bounding(Eigen::MatrixXd& input);
	
	void RGradient(Eigen::MatrixXd& Rgrad);	
	void X_vGradient(Eigen::MatrixXd& Xgrad, int v);
	void W_bGradient(Eigen::MatrixXd& Wbgrad, Eigen::MatrixXd& Constbgrad,int b);

	Eigen::MatrixXi G;	// adjacency matrix
	Eigen::MatrixXi A; // attribute matrix
	// latent variable. See the paper to know each meaning
	Eigen::MatrixXd X;
	Eigen::MatrixXd R;
	Eigen::MatrixXd W;
	Eigen::MatrixXd Constant; // this term is used to store regression constant
	Eigen::MatrixXd x_sum_vector;

	int k_target; // number of roles
	double alpha_r, alpha_x, alpha_w; // coefficient of regularization. See the paper
	double beta_x, beta_r, beta_w; //used for zoom. The degree of closeness enforcement
	double alpha; // weight on the attribute term in cost function
	bool digraph; // determine if it is directed graph
};

#endif