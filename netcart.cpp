#include "netcart.h"
#include "util.h"


Netcart::Netcart()
{
	k_target = 0;	
	alpha_r = 0;
	alpha_x = 0;
	alpha_w = 0;
	alpha = 0;
	beta_x = 0;
	beta_r = 0;
	beta_w = 0;
}

Netcart::Netcart(
	string graph_file, string attri_file,  int _k_target, 
	double _alpha, bool _digraph)
{
	SetNumRole(_k_target);
	SetAttriWeight(_alpha);
	SetDigraph(_digraph);
	SetAttributedGraph(graph_file, attri_file);
}

void 
Netcart::SetRegularization(double _alpha_x, double _alpha_r, double _alpha_w)
{
	alpha_x = _alpha_x;
	alpha_w = _alpha_w;
	alpha_r = _alpha_r;
}

void 
Netcart::SetLearningRate(double _beta_x, double _beta_r, double _beta_w)
{
	beta_x = _beta_x;
	beta_r = _beta_r;
	beta_w = _beta_w;
}

void 
Netcart::SetAttributedGraph(string graph_file, string attri_file)
{
	vector< pair<int,int> > EdgeList;
	set<int> NodeList;
	vector< pair<int,int> > NodeAttribute;
	set<int> AttributeList;	

	//reading in graph file
	fstream myfile;
	myfile.open(graph_file.c_str(),fstream::in);
	while(!myfile.eof())
	{
		int node1, node2;
		myfile >> node1 >> node2;
		EdgeList.push_back(make_pair(node1, node2));
		NodeList.insert(node1);
		NodeList.insert(node2);
	}
	myfile.close();

	size_t NumOfNode = NodeList.size();
	G = Eigen::MatrixXi::Zero(NumOfNode,NumOfNode);
	for (auto edge : EdgeList)
		G(edge.first, edge.second) = 1;

	//reading in attribute file
	myfile.open(attri_file.c_str(), fstream::in);
	while(!myfile.eof())
	{
		int node, attribute;
		myfile >> node >> attribute;
		NodeAttribute.push_back(make_pair(node, attribute));
		AttributeList.insert(attribute);
	}

	size_t NumOfAttri = AttributeList.size();
	A = Eigen::MatrixXi::Zero(NumOfNode, NumOfAttri);
	for (auto has_attri : NodeAttribute)
		A(has_attri.first, has_attri.second) = 1;

	cout << "There are " << NumOfNode << " nodes, " << NumOfAttri 
	<< " attributes each node and " << EdgeList.size() << " edges." << endl;
	//set the size of X and W
	X = Eigen::MatrixXd::Zero(k_target, NumOfNode);
	W = Eigen::MatrixXd::Zero(k_target, NumOfAttri);
	Constant = Eigen::MatrixXd::Zero(1, NumOfAttri);
	R = Eigen::MatrixXd::Zero(k_target, k_target);
}

void Netcart::SetNumRole(int _k_target ){k_target = _k_target;}
void Netcart::SetAttriWeight(double _alpha){alpha = _alpha;}
void Netcart::SetDigraph(bool _digraph){digraph = _digraph;}

void Netcart::initialization()
{
	srand(13);
	R.setRandom();
	R = R + Eigen::MatrixXd::Ones(R.rows(),R.cols());//non-negative constraint
	if(!digraph)
	{
		Eigen::MatrixXd M = R; 
		R = (M + M.transpose())/2;
	}
	X.setRandom();
	X = X + Eigen::MatrixXd::Ones(X.rows(),X.cols());//non-negative constraint
	x_sum_vector = X.rowwise().sum();
	W.setRandom();
	Constant.setRandom();
}

double Netcart::LogLikelihood()
{
	return LogLikelihoodGraph() + LogLikelihoodAttri();
}

double Netcart::LogLikelihoodGraph()
{
	double result = 0;
	for (int u = 0; u < G.rows(); ++u)
	{
		if(!digraph)
		{
			for (int v = 0; v < G.cols(); ++v)
			{
				if(u == v) continue;
				double predict = predictEdge(u,v);
				if (G(u,v) == 0)
				{
					if(predict == 1)
						predict = 1 - SMALLEST_NUM;
					result += log(1-predict);
				}
				else
				{
					if(predict == 0)
						predict = SMALLEST_NUM;
					result += log(predict);
				}
			}
		}
		else
		{
			for (int v = u+1; v < G.cols(); ++v)
			{
				double predict = predictEdge(u,v);
				if (G(u,v)==0)
					result += log(1-predict);
				else
					result += log(predict);
			}
		}	
	}
	return result;
}

double Netcart::LogLikelihoodAttri()
{
	double result = 0;
	for (int v = 0; v < A.rows(); ++v)
	{
		for (int i = 0; i < A.cols(); ++i)
		{
			double predict = predictAttri(v,i);
			if (A(v,i)==0)
			{
				if (predict == 1)
					predict = 1 - SMALLEST_NUM;
				result += log(1-predict);
			}
			else
			{
				if (predict == 0)
					predict = SMALLEST_NUM;
				result += log(predict);
			}
		}
	}
	return result;
}

double Netcart::predictEdge(int u, int v)
{
	//calculating the rho_{uv}
	double rho = X.col(u).transpose() * R* X.col(v);
	if(rho < 0)
	{
		cout << "predict edge probability less than 0!" << endl; 
		exit(0);
	}
	return 1 - exp(-rho);
}

double Netcart::predictAttri(int v, int i)
{
	double miu = W.col(i).transpose()*X.col(v) + Constant(0, i);
	return 1/(1+exp(-miu));
}

void Netcart::PrintXRW()
{
	PrintX();
	PrintR();
	PrintW();
}
void Netcart::PrintX()
{
	cout << "X: " << endl;
	cout << X << endl;
}
void Netcart::PrintR()
{
	cout << "R: " << endl;
	cout << R << endl;
}
void Netcart::PrintW()
{
	cout << "W: " << endl;
	cout << W << endl;
}

void Netcart::SaveXRW(string path)
{
	SaveX(path);
	SaveR(path);
	SaveW(path);
	SaveConst(path);
}

void Netcart::SaveConst(string path){WriteMat(Constant, path+"Const.txt");}
void Netcart::SaveX(string path){WriteMat(X, path+"X.txt");}
void Netcart::SaveR(string path){WriteMat(R, path+"R.txt");}
void Netcart::SaveW(string path){WriteMat(W, path+"W.txt");}

double Netcart::CostFunction()
{
	return (1-alpha)*LogLikelihoodGraph() + alpha*LogLikelihoodAttri()
			- alpha_r*L1Norm(R) - alpha_x*L1Norm(X) - 
			alpha_w*L1Norm(W) - alpha_w*L1Norm(Constant);
}

void Netcart::Optimize(int max_iterate_x, int max_iterate_r, int max_iterate_w)
{
	double CostCurrent = CostFunction();
	cout << "starting at: " << CostCurrent << endl;
	for (int i = 0; i < MAX_ITERATE; ++i)
	{
		OptimizeR(max_iterate_r);
		OptimizeX(max_iterate_x);
		OptimizeW(max_iterate_w);
		double CostOld = CostCurrent;
		CostCurrent = CostFunction();
		cout << "at " << i+1 << "th iteration: " << CostCurrent << endl;
		if (Converge(CostOld, CostCurrent))
			break;
	}
	cout << "saving results..." << endl;
	SaveXRW();
	return 0;
}

// REQUIRE: R/X/Wgrad has to be zero matrix
void Netcart::RGradient(Eigen::MatrixXd& Rgrad)
{
	Eigen::MatrixXd res_edge = Rgrad; 
	Eigen::MatrixXd XXsumEdge = Rgrad;
	for (int v = 0; v < G.rows(); ++v)
	{
		XXsumEdge += X.col(v)*X.col(v).transpose();
		for (int u = 0; u < G.cols(); ++u)
		{
			if (u == v)
			{
				continue;
			}
			if(G(v,u) == 1)
			{
				double temp = predictEdge(v,u);
				if (temp == 0)
				{
					temp = SMALLEST_NUM;
				}
				if (temp == 1)
				{
					temp = 1 - SMALLEST_NUM;
				}
				Eigen::MatrixXd tmp = X.col(v)*X.col(u).transpose();
				res_edge += tmp*(1-temp)/temp;
				XXsumEdge += tmp;
			}
		}
	}
	Eigen::MatrixXd res_non_edge =  x_sum_vector * x_sum_vector.transpose() - XXsumEdge;
	Rgrad = (1 - alpha)*(res_edge - res_non_edge) - alpha_r * Sign(R);
	normalization(Rgrad);
}


void Netcart::X_vGradient(Eigen::MatrixXd& Xgrad, int v)
{
	Eigen::MatrixXd Xvgrad_edge	= Eigen::MatrixXd::Zero(X.rows(),1);
	if (digraph)
	{
		Eigen::MatrixXd	out_nbr_v = Eigen::MatrixXd::Zero(X.rows(),1);
		Eigen::MatrixXd out_Xsum_v = Eigen::MatrixXd::Zero(X.rows(),1);
		Eigen::MatrixXd	in_nbr_v = Eigen::MatrixXd::Zero(X.rows(),1);
		Eigen::MatrixXd in_Xsum_v = Eigen::MatrixXd::Zero(X.rows(),1);
		for (int u = 0; u < G.cols(); ++u)
		{
			if (v == u)
			{
				continue;
			}
			// out-neighbour of v
			if(G(v,u) == 1)
			{
				double temp = predictEdge(v,u);
				if (temp == 0)
				{
					temp = SMALLEST_NUM;
				}
				if (temp == 1)
				{
					temp = 1 - SMALLEST_NUM;
				}
				out_nbr_v +=  (1 - temp) * X.col(u)/temp;
				out_Xsum_v += X.col(u);
			}
			// in-neighbour of v
			if (G(u,v) == 1)
			{
				double temp = predictEdge(u,v);
				cout << u << ": " << temp << endl;
				if (temp == 0)
				{
					temp = SMALLEST_NUM;
				}
				if (temp == 1)
				{
					temp = 1 - SMALLEST_NUM;
				}
				in_nbr_v +=  (1 - temp) * X.col(u)/temp;
				in_Xsum_v += X.col(u);
			}
		}
		Eigen::MatrixXd out_non_nbr_v = x_sum_vector - out_Xsum_v - X.col(v);
		Eigen::MatrixXd in_non_nbr_v = x_sum_vector - in_Xsum_v - X.col(v);		
		Xvgrad_edge = R*(out_nbr_v - out_non_nbr_v) + R.transpose()*(in_nbr_v - in_non_nbr_v);
	}
	Eigen::MatrixXd	Xvgrad_attri = Eigen::MatrixXd::Zero(X.rows(),1);
	for (int i = 0; i < A.cols(); ++i)
	{
		double temp = predictAttri(v,i); //sigma
		if (A(v,i) == 1)
		{
			if (temp == 1)
			{
				temp = 1 - SMALLEST_NUM;
			}
			Xvgrad_attri += (1 - temp) * W.col(i);
		}
		else
		{
			if (temp == 0)
			{
				temp = SMALLEST_NUM;
			}
			Xvgrad_attri -= temp * W.col(i);
		}
	}
	//add regularization of Xgrad
	Eigen::MatrixXd X_v = X.col(v);
	Eigen::MatrixXd Xvgrad = (1 - alpha) * Xvgrad_edge + alpha * Xvgrad_attri - alpha_x * Sign(X_v);
	normalization(Xvgrad);
	Xgrad.col(v) = Xvgrad;
}

void Netcart::W_bGradient(Eigen::MatrixXd& Wgrad, Eigen::MatrixXd &Constgrad, int b)
{
	Eigen::MatrixXd W_b = W.col(b);
	Eigen::MatrixXd Const_b = Constant.col(b);
	Eigen::MatrixXd Wbgrad = Eigen::MatrixXd::Zero(Wgrad.rows(), 1);
	Eigen::MatrixXd Constbgrad = Eigen::MatrixXd::Zero(1, 1);
	// just a work around for the "+=" overload of Eigen::Dense.
	// Eigen is terrible for developer
	Eigen::MatrixXd one = Eigen::MatrixXd::Ones(1, 1);
	Eigen::MatrixXd attrib_term = Eigen::MatrixXd::Zero(Wgrad.rows(), 1);
	for (int v = 0; v < A.rows(); ++v)
	{
		double temp = predictAttri(v,b); 
		if (temp == 0)
		{
			temp = SMALLEST_NUM;
		}
		if (temp == 1)
		{
			temp = 1 - SMALLEST_NUM;
		}
		if (A(v,b) == 1)
		{
			attrib_term += ( 1 - temp ) * X.col(v);
			Wbgrad += (1 - temp) * X.col(v);
			Constbgrad += (1 - temp) * one;
		}
		else
		{
			Wbgrad -= temp * X.col(v);
			Constbgrad -= temp * one;
		}
	}
	Wbgrad *= alpha;
	Constbgrad *= alpha;
	// regularization
	Wbgrad -= alpha_w * Sign(W_b);
	Constbgrad -= alpha_w * Sign(Const_b);
	if (alpha != 0)
	{
		normalization(Wbgrad);
	}
	Wgrad.col(b) = Wbgrad;
	Constgrad.col(b) = Constbgrad;
}

bool Netcart::Converge(double CostOld, double CostCurrent)
{
	return abs((CostOld - CostCurrent)/CostOld) <= 1e-5;
}


void
Netcart::Bounding(Eigen::MatrixXd &input)
{
	for (int i = 0; i < input.rows(); ++i)
	{
		for (int j = 0; j < input.cols(); ++j)
		{
			if (input(i,j) < 0)
			{
				input(i,j) = SMALLEST_NUM;
			}
			if (input(i,j) > 7)
			{
				input(i,j) = 7 - SMALLEST_NUM;
			}
		}
	}
}