#include "netcart.h"
#include "util.h"

// TODO: the converging speed is not good. 
// TODO: add multi-threading or exploiting Eigen.
// TODO: add zoom

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
//	srand(time(NULL));
	R.setRandom();
	R = R + Eigen::MatrixXd::Ones(R.rows(),R.cols());//non-negative constraint
	if(!digraph)
	{
		Eigen::MatrixXd M = R; 
		R = (M + M.transpose())/2;
	}
	X.setRandom();
	X = X + Eigen::MatrixXd::Ones(X.rows(),X.cols());//non-negative constraint
	W.setRandom();
	Constant.setRandom();
}

double Netcart::LogLikelihood()
{
	return LogLikelihoodGraph() + LogLikelihoodAttri();
}

double Netcart::LogLikelihoodGraph()
{
	double resedge = 0;
	x_sum_vector = X.rowwise().sum();
	double nonresedge = (x_sum_vector.transpose() * R * x_sum_vector)(0);
	for (int u = 0; u < G.rows(); ++u)
	{
		nonresedge -= X.col(u).transpose() * R * X.col(u);
	}
	if(!digraph)
	{
		nonresedge /= 2;
	}
	for (int v = 0; v < G.cols(); ++v)
	{
		for (int u = 0; u < G.cols(); ++u)
		{
			if(u == v) continue;
			if (G(u,v) == 1)
			{
				double predict = X.col(u).transpose() * R * X.col(v);
				if(predict == 0)
					predict = ZERO;
				nonresedge -= predict;
				resedge += log(1 - exp(-predict));
			}
		}
	}
	nonresedge *= -1;
//	cout << "resedge: " << resedge << " nonresedge: " << nonresedge << endl;
//	cout << "edge: " << resedge + nonresedge << endl;
	return resedge + nonresedge;
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
					predict = 1 - ZERO;
				result += log(1-predict);
			}
			else
			{
				if (predict == 0)
					predict = ZERO;
				result += log(predict);
			}
		}
	}
//	cout << "attribute: " << result << endl;
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

void Netcart::Optimize(int max_iterate_x, int max_iterate_r, int max_iterate_w, string path)
{
	double CostCurrent = CostFunction();
	cout << "starting at: " << CostCurrent << endl;
	for (int i = 0; i < MAX_ITERATE; ++i)
	{
		OptimizeW(max_iterate_w);
		OptimizeR(max_iterate_r);
		if (!digraph)
		{
			R = R + R.transpose()/2;
		}
		OptimizeX(max_iterate_x);

		double CostOld = CostCurrent;
		CostCurrent = CostFunction();
		cout << "at " << i+1 << "th iteration: " << CostCurrent << endl;
		if (Converge(CostOld, CostCurrent))
			break;
	}
	cout << "saving results..." << endl;
	SaveXRW(path);
}

void Netcart::OptimizeR(int max_iterate_r)
{
	x_sum_vector = X.rowwise().sum();
	double CostCurrent = CostFunction();
	for (int i = 0; i < max_iterate_r; ++i)
	{
		Eigen::MatrixXd Rgrad = Eigen::MatrixXd::Zero(R.rows(), R.cols());
		double CostOld = CostCurrent;
		RGradient(Rgrad);
		R += beta_r * Rgrad;
		Bounding(R);
		CostCurrent = CostFunction();
		if (Converge(CostOld, CostCurrent))
			break;	
	}
}

void Netcart::OptimizeX(int max_iterate_x)
{
	x_sum_vector = X.rowwise().sum();	
	double CostCurrent = CostFunction();
	for (int v = 0; v < G.cols(); ++v)
	{
		for (int i = 0; i < max_iterate_x; ++i)
		{
			x_sum_vector = X.rowwise().sum();	
			Eigen::MatrixXd X_vgrad = Eigen::MatrixXd::Zero(X.rows(), 1);
			Eigen::MatrixXd X_v = X.col(v);
			double CostOld = CostCurrent;
			X_vGradient(X_vgrad, v);
			X_v += beta_x * X_vgrad;
			Bounding(X_v);
			X.col(v) = X_v;
			CostCurrent = CostFunction();
			if (Converge(CostOld, CostCurrent))
				break;
		}
	}
}

void Netcart::OptimizeW(int max_iterate_w)
{
	double CostCurrent = CostFunction();
	for (int b = 0; b < W.cols(); ++b)
	{
		for (int i = 0; i < max_iterate_w; i++)
		{
			Eigen::MatrixXd W_b = W.col(b);
			Eigen::MatrixXd Const_b = Constant.col(b);
			Eigen::MatrixXd W_bgrad = Eigen::MatrixXd::Zero(W.rows(), 1);
			Eigen::MatrixXd Const_bgrad = Eigen::MatrixXd::Zero(Constant.rows(), 1);
			double CostOld = CostCurrent;
			W_bGradient(W_bgrad, Const_bgrad, b);
			W_b += beta_w * W_bgrad;
			W.col(b) = W_b;
			Const_b	+= beta_w * Const_bgrad;
			Constant.col(b) = Const_b;
			CostCurrent = CostFunction();
			if (Converge(CostOld, CostCurrent))
				break;
		}
	}
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
					temp = ZERO;
				}
				if (temp == 1)
				{
					temp = 1 - ZERO;
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


void Netcart::X_vGradient(Eigen::MatrixXd& X_vgrad, int v)
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
					temp = ZERO;
				}
				if (temp == 1)
				{
					temp = 1 - ZERO;
				}
				out_nbr_v +=  (1 - temp) * X.col(u)/temp;
				out_Xsum_v += X.col(u);
                    
                //in_nbr_v += (1 - temp) * X.col(v)/temp;
                //in_Xsum_v += X.col(v);
			}
			// in-neighbour of v
	    	if (G(u,v) == 1)
	    	{
	    		double temp = predictEdge(u,v);
	    		if (temp == 0)
	    		{
	    			temp = ZERO;
	    		}
	    		if (temp == 1)
	    		{
	    			temp = 1 - ZERO;
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
				temp = 1 - ZERO;
			}
			Xvgrad_attri += (1 - temp) * W.col(i);
		}
		else
		{
			if (temp == 0)
			{
				temp = ZERO;
			}
			Xvgrad_attri -= temp * W.col(i);
		}
	}
	//add regularization of Xgrad
	Eigen::MatrixXd X_v = X.col(v);
	X_vgrad = (1 - alpha) * Xvgrad_edge + alpha * Xvgrad_attri - alpha_x * Sign(X_v);
	normalization(X_vgrad);
}

void Netcart::W_bGradient(Eigen::MatrixXd& W_bgrad, Eigen::MatrixXd &Const_bgrad, int b)
{
	Eigen::MatrixXd W_b = W.col(b);
	Eigen::MatrixXd Const_b = Constant.col(b);
	// just a work around for the "+=" overload of Eigen::Dense.
	// Eigen is terrible for developer
	Eigen::MatrixXd one = Eigen::MatrixXd::Ones(1, 1);
	for (int v = 0; v < A.rows(); ++v)
	{
		double temp = predictAttri(v,b); 
		if (temp == 0)
		{
			temp = ZERO;
		}
		if (temp == 1)
		{
			temp = 1 - ZERO;
		}
		if (A(v,b) == 1)
		{
			W_bgrad += (1 - temp) * X.col(v);
			Const_bgrad += (1 - temp) * one;
		}
		else
		{
			W_bgrad -= temp * X.col(v);
			Const_bgrad -= temp * one;
		}
	}
	W_bgrad *= alpha;
	Const_bgrad *= alpha;
	// regularization
	W_bgrad -= 2 * alpha_w * W_b;
	Const_bgrad -= 2 * alpha_w * Const_b;
	if (alpha != 0)
	{
		normalization(W_bgrad);
	}
}

bool Netcart::Converge(double CostOld, double CostCurrent)
{
	return (CostCurrent - CostOld)/abs(CostOld) <= 1e-5;
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
				input(i,j) = ZERO;
			}
			if (input(i,j) > 7)
			{
				input(i,j) = 7 - ZERO;
			}
		}
	}
}
