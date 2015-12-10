/*-------------------------------------------------------------
 * Author: Yuchen Lu
 * Compilation: make
 * Execution: the program is under bin/
 * 
 * Perform a RoDan algorithm on hospital data set with k equals 
 * to 4, and save the result to current directory.
 *
 *-----------------------------------------------------------*/

#ifndef NETCART_H
#define NETCART_H

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

const int MAX_ITERATE = 1;
const double ZERO = 1e-14;

using namespace std;

class Netcart
{
public:
    /**
     * default constructor.
     */
    Netcart();

    /**
     * graph_file: the location of .edge file.
     * attri_file: the location of .nodefeat file.
     * _k_target: the number of roles to be considered.
     * _alpha: the weight of attribute term in cost function. Between 0 to 1.
     * _digraph: True if the graph is directed. False otherwise.
     */
	Netcart(string graph_file, string attri_file, int _k_target = 4, 
			double _alpha = 0.5, bool _digraph = true);

    /**
     * The parameters for L1 norm regularization.
     */
	void SetRegularization(double alpha_x = 0.1, double alpha_r = 0.9375, double alpha_w = 0.0001875);

    /**
     * A bunch of functions used for set data member.
     */
	void SetAttributedGraph(string graph_file, string attri_file);
    void SetNumRole(int _k_target);
    void SetAttriWeight(double _alpha);
	void SetDigraph(bool _digraph);
	void SetLearningRate(double _beta_x, double _beta_r, double _beta_w);
	void SetX(Eigen::MatrixXd _X){X = _X;}
	void SetR(Eigen::MatrixXd _R){R = _R;}
	void SetW(Eigen::MatrixXd _W){W = _W;}

    /**
     * A bunch of functions to get the data member.
     */
	void PrintXRW();
	void PrintX();
	void PrintR();
	void PrintW();
	void SaveXRW(string path = "");
	void SaveX(string path = "");
	void SaveR(string path = "");
	void SaveW(string path = "");
	void SaveConst(string path = "");

    /**
     * Initialize X, R, and W. Should be called before optimization.
     */
	void initialization();

    /**
     * likelihood computation
     */
	double LogLikelihood();
	double LogLikelihoodGraph();
	double LogLikelihoodAttri();

    /**
     * The cost function includes total likelihood and regulariation
     */
	double CostFunction();
    
    /**
     * max_iterate_x/r/w: the number of iteration for X/R/W. We use a group coordinates
     * gradient ascend here.
     * path: the location to save the result.
     */
	void Optimize(int max_iterate_x, int max_iterate_r, int max_iterate_w, string path = "");
	void OptimizeR(int max_iterate_r);
	void OptimizeX(int max_iterate_x);
	void OptimizeW(int max_iterate_w);

private:
    /**
     * Returns the probabilty that there is an edge (i,j).
     */
	double predictEdge(int i, int j);

    /**
     * Returns the probability that node v has attribute i.
     */
	double predictAttri(int v, int i);
	
    /**
     * Converging condition to stop the iteration.
     */
    bool Converge(double oldCost, double newCost);

    /**
     * This used to bound the range of our X/R/W during optimization
     */
	void Bounding(Eigen::MatrixXd& input);
	
    /**
     * Returns the gradient for R
     */
	void RGradient(Eigen::MatrixXd& Rgrad);	
	
    /**
     * Returns the gradient for X.col(v). 
     * Notice: Xgrad still has the same size of X, but
     * only column v is changed.
     */
    void X_vGradient(Eigen::MatrixXd& Xgrad, int v);

    /**
     * Returns the gradient of W.col(b)
     * The same notice as X_vGradient.
     */
	void W_bGradient(Eigen::MatrixXd& Wbgrad, Eigen::MatrixXd& Constbgrad,int b);

    /**
     * Adjacency Matrix. 
     * G(i,j) means there is edge from i to j.
     */
	Eigen::MatrixXi G;

    /**
     * Attribute Matrix.
     * A(i,j) means node i has attribute j
     */
	Eigen::MatrixXi A;

	/**
     * latent parameters to be optimized. 
     * See the paper for detailed explanation.
     */
	Eigen::MatrixXd X;
	Eigen::MatrixXd R;
	Eigen::MatrixXd W;

    /**
     * This is for logistic regression constant.
     */
	Eigen::MatrixXd Constant;

    /**
     * This term is cached before each iteration.
     */
	Eigen::MatrixXd x_sum_vector;

    /**
     * Number of Roles.
     */
	int k_target; 

    /**
     * Regularization coefficients.
     */
	double alpha_r, alpha_x, alpha_w; 

    /**
     * Learning rate
     */
    double beta_x, beta_r, beta_w;

    /**
     * Weight attribute.
     */
	double alpha;

    /**
     * Directed graph indicator
     */
	bool digraph;
};

#endif
