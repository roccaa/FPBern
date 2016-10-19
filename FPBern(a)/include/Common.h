/*
 * Common.hpp
 *
 *  Created on: Jul 25, 2014
 *      Author: dreossi
 */

#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <ostream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <fstream>

#include <math.h>
#include <cmath>

#include <vector>

#include "assert.h"


#include <ginac/ginac.h>
#include <boost/lexical_cast.hpp>

//#define epsilon 1e-30

using namespace std;
using namespace GiNaC;

typedef std::vector<double> doubleVector;

typedef std::vector<int> intVector;

typedef vector < vector <double> > vector_matrix;

typedef int node_type;
enum {IS_NODE, IS_LEAF};

struct synthesizer_opt{
	bool largest_para_set;	// largest parameter set
};

struct sapo_option{
	int trans;				// transformation (0: static, 1: dynamic)
	double alpha;			// decomposition weight
	int decomp;				// number of decompositions (0: none, >0: yes)
	string plot;			// the name of the file were to plot the reach set
	bool verbose;			// display info
};

struct approx_infos{
	string tree_ID;				// Tree containning the partition
	string leave_ID;			// partition where we add the approx
	vector< std::pair<std::string,ex> > list_approx; // list of approx to add: a symbolic variable, and its local value;
	vector< ex > list_approx_ex; // list of approx in the form 'approx == expression' as ginac expresssion;
};


struct poly_values{
	vector<double> base_vertex;
	vector< vector<double> > versors;
	vector<double> lenghts;
};

typedef int equation_type;
enum {Polynomial_ONLY,RationalPolynomial_ONLY,MixedEquation};

typedef int bernstein_computation_method;
enum {IMPLICIT_SYM,IMPLICIT_NUM,EXPLICIT_SYM};

struct bern_info
{
	bernstein_computation_method method;
	// Si la liste des coefficients est partielle ou non
	bool partial;
//	map<intVector,lst> sym_coeff_map;
	vector<ex>* coeff_numer;
	vector<ex>* coeff_denom;

	vector<ex>* min_coeffs_guess;
	vector<ex>* max_coeffs_guess;

	vector<double>* coeff_min;
	vector<double>* coeff_max;
	double min;
	double max;
};

enum {OFO,AFO};

enum {PARTIAL_TEMPLATE_GUESS,FULL_TEMPLATE_GUESS};

typedef int representation;
enum {BOX,PARALLELOTOPE,BUNDLE};

struct bernstein_options
{
	bool rationnal;
	bernstein_computation_method method;
};

struct reach_info
{
	int dim;
	vector_matrix L;
	vector<intVector> T;
	representation type;
	doubleVector offmin;
	doubleVector offmax;
};

enum formula_type {ATOM,CONJUNCTION,DISJUNCTION,UNTIL,ALWAYS,EVENTUALLY,TEMPO_ATOM};

enum {DECREASING,INCREASING};

typedef int numerical_algorithm;
enum {EULER,RUNGE_KUTTA};

typedef int splitting_decision_method;
enum {Full_Hessian_Full_Box, Partial_Hessian_Box_GridBased};

typedef int splitting_direction_method;
enum {MaxJacobian_FULL_LP,MaxJacobian_FULL_NLP,Max_Box_Size};


struct reach_options
{
	splitting_decision_method dec_meth;
	splitting_direction_method dir_meth;
	bernstein_computation_method compu_meth;
	int number_of_steps;
	numerical_algorithm sim_algo;
	int sim_sampling;
	double tstep;
	int max_split;
	string time_unit;
};

typedef int splitting_type;
enum {no_splitting,virtual_splitting,real_splitting};

enum {one_splitting,partial_splitting,full_splitting};

enum {C_FUNCTION,GINAC_EX};

enum {STATIC,DYNAMIC};



inline std::ostream& operator<<(std::ostream& out, doubleVector v) {
	cout.precision(8);
	for (unsigned int i = 0; i < v.size(); i++)
		out << v[i] << " ";
	out << std::endl;
	return out;
}

inline std::ostream& operator<<(std::ostream& out, intVector v) {
	for (unsigned int i = 0; i < v.size(); i++)
		out << v[i] << " ";
	out << std::endl;
	return out;
}

inline std::ostream& operator<<(std::ostream& out, vector_matrix v) {
	for(unsigned int i = 0; i < v.size(); i++)
	{
		for(unsigned int j = 0; j < v[i].size(); j++)
		{
			out << v[i][j] << " ";
		}
	out << std::endl;
	}
	cout << "done!\n";
	return out;
}

inline std::ostream& operator<<(std::ostream& out, vector< vector<int> > v) {
	for(unsigned int i = 0; i < v.size(); i++)
	{
		for(unsigned int j = 0; j < v[i].size(); j++)
		{
			out << v[i][j] << " ";
		}
	out << std::endl;
	}
	cout << "done!\n";
	return out;
}

inline  lst vect_to_lst(vector<ex> vec)
{
	lst res;
	for(int i=0;i<vec.size();i++)
	{
		res.append(vec[i]);
	}
	return res;
}

inline intVector index_to_mindex(int index, intVector shift_table)
{
	intVector res = intVector(shift_table.size(),0);
	doubleVector table = doubleVector(shift_table.size(),1);
	double total = shift_table[0];
	for(int i=table.size()-2;i>=0;i--)
	{
		table[i] = total;
		total = total*shift_table[i];
	}
	int dim = shift_table.size();
	int retenu = index;
	for(int i=0;i<dim;i++){
		//cout << "###retenu = " << retenu << endl;
		//cout << "shift = " << ((int)table[i]) << endl;
		if(i==dim-1)
			res[i] = retenu;
		else{
			res[i] = retenu/((int)table[i]);
			retenu = retenu%((int)table[i]);
		}
		//cout << "index[i] = " << res[i] << endl;
	}

	return res;

}

//inline arma::Mat<double> normalize_row(arma::Mat<double> m)
//{
//	arma::Mat<double> res = m;
//	cout << "in normalize_row:\n";
//	for(int i=0;i<m.n_rows;i++)
//	{
//		cout << "norm = " << norm(res.row(i),2) <<"\n";
//		res.row(i) = res.row(i)/norm(res.row(i),2);
//	}
//	//res.print();
//	return res;
//}

inline double **mat_create(int rows, int cols)
{
	double **M; int i;

	M = (double**) malloc(rows * sizeof(double*));
	for (i=0; i<rows; i++)
		M[i] = (double*) malloc(cols * sizeof(double));

	return M;
}

//inline double ** arma_to_Cmat(arma::Mat<double>* arma_m)
//{
//
//	int rows = arma_m->n_cols;
//	int cols = arma_m->n_rows;
//	double **M;
//
//	M = (double**) malloc(rows * sizeof(double*));
//	for (int i=0; i<rows; i++)
//		M[i] = (double*) malloc(cols * sizeof(double));
//
//	for(int i=0; i<rows; i++)
//	{
//		for (int j=0; j<cols; j++)
//		{
//			M[i][j] = arma_m->at(j,i);
//		}
//	}
//
//	return M;
//}


//inline vector_matrix arma2std_Mat(arma::mat arma_m)
//{
//	assert(arma_m.n_elem != 0);
//	vector_matrix res;
//	for(int i=0;i<arma_m.n_rows;i++)
//	{
//		doubleVector r = doubleVector(arma_m.n_cols);
//		for(int j=0;j<arma_m.n_cols;j++)
//		{
//			r[j] = arma_m(i,j);
//		}
//		res.push_back(r);
//	}
//	return res;
//}

//inline arma::Mat<double> std2arma_Mat(vector_matrix vect_m)
//{
//	assert(vect_m.size() != 0);
//	arma::Mat<double> res = arma::Mat<double>(vect_m.size(),vect_m[0].size());
//	for(int i=0;i<vect_m.size();i++)
//	{
//		arma::Row<double> r = arma::Row<double>(vect_m[i]);
//		res.row(i) = r;
//
//	}
//	return res;
//}


inline bool operator==(doubleVector a, doubleVector b)
{
	if(a.size() != b.size())
		return false;
	for(int i=0;i<a.size();i++)
	{
		if(a[i]!=b[i])
			return false;
	}
	return true;
}


inline int intersect_matrix(vector_matrix *A,doubleVector *b,doubleVector constraint,double offset)
{
	for(int i=0;i<A->size();i++)
	{
		if(A->at(i) == constraint)
		{
			b->at(i) = min(b->at(i),offset);
			return i;
		}
	}
	A->push_back(constraint);
	b->push_back(offset);
	return -2;
}

inline double read_double(string *input,string sep)
{
	int i = input->find(sep);
	string res = input->substr(0,i);
	input->erase(input->begin(),input->begin()+i+1);
	double res_d;
	try
	{
		//cout << "string = " << res << endl;
		res_d = boost::lexical_cast<double>(res);
	}
	catch (boost::bad_lexical_cast const&)
	{
		cout << "Wrong input - not a double\n";
	   	assert(0);
	}
	return res_d;
}




inline std::istream& operator>>(std::istream& str, std::string line)
{
	std::getline(str,line);
    return str;
}

inline doubleVector read_doubleVector(string *input,string sep,int d)
{
	doubleVector res = doubleVector(d);
	//cout << *(input) << endl;
	for(int i=0;i<d;i++)
	{
		double v = read_double(input,sep);
		res[i] = v;
	}
	return res;
}

inline string read_string(string *input,string sep)
{
	int i = input->find(sep);
	string res = input->substr(0,i);
	input->erase(input->begin(),input->begin()+i+1);
	return res;
}

//inline double compute_angle(doubleVector d1,doubleVector d2)
//{
//	arma::vec a = arma::vec(d1);
//	arma::vec b = arma::vec(d2);
//	double s = arma::norm_dot(a,b);
//
//	return acos(s);
//}

//inline vector<double> compute_listAngles(vector<doubleVector> list_directions)
//{
//	doubleVector list_angles;
//	for(int i=0;i<list_directions.size();i++)
//	{
//		for(int j=0;j<list_directions.size();j++)
//		{
//			double alpha = compute_angle(list_directions[i],list_directions[j]);
//			list_angles.push_back(alpha);
//		}
//	}
//
//	return list_angles;
//}

inline vector<double> * copy(vector<double> * a){
	vector<double> * b = new vector<double>(a->size(),0);
	for(int i=0;i<a->size();i++){
		(*b)[i] = (*a)[i];
	}
	return b;
}

inline int max(const vector<int>& v) {
	assert(v.size() > 0);
	int max = v[0];
	for (unsigned i = 1; i < v.size(); i++)
		if (max < v[i])
			max = v[i];
	return max;
}

inline double max(const vector<double>& v) {
	assert(v.size() > 0);
	double max = v[0];
	for (unsigned i = 1; i < v.size(); i++)
		if (max < v[i])
			max = v[i];
	return max;
}

inline doubleVector operator+(doubleVector v1, doubleVector v2){
	doubleVector res = doubleVector(v1.size());
	assert(v1.size() == v2.size());
	for(int i=0;i<v1.size();i++){
		res[i] = v1[i]+v2[i];
	}
	return res;
}

inline doubleVector operator-(doubleVector v1, doubleVector v2){
	doubleVector res = doubleVector(v1.size());
	assert(v1.size() == v2.size());
	for(int i=0;i<v1.size();i++){
		res[i] = v1[i]-v2[i];
	}
	return res;
}
inline doubleVector operator*(doubleVector v1, double var){
	doubleVector res = doubleVector(v1.size());
	for(int i=0;i<v1.size();i++){
		res[i] = v1[i]*var;
	}
	return res;
}
inline doubleVector operator*(double var, doubleVector v1){
	doubleVector res = doubleVector(v1.size());
//	cout << "var = " << var <<endl;
//	cout << "v1 = " << v1;
	for(int i=0;i<v1.size();i++){
		res[i] = v1[i]*var;
	}
//	cout << "res = " << res;
	return res;
}


#endif /* COMMON_HPP_ */
