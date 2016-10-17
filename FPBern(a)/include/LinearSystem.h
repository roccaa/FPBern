/*
 * LinearSystem.h
 *
 *  Created on: Oct 24, 2014
 *      Author: Tommaso Dreossi
 */

#ifndef LINEARSYSTEM_H_
#define LINEARSYSTEM_H_

#include "Common.h"
#include <glpk.h>

class LinearSystem {

private:
	int n_vars;			// number of variables
	lst vars;			// 	list of variables
	lst constraints;	//	list of constraints
	vector< vector<double> > A; // matrix A
	vector< double > b; 		// vector b

	bool isIn(vector< double > Ai, double bi);	// check if a constraint is already in
	void initLS();				// initialize A and b
	double solveLinearSystem(vector< vector< double > > A, vector< double > b, vector< double > obj_fun, int min_max);
	bool zeroLine(vector<double> line);


public:

	// For splitting use only
	vector_matrix split_dir;// trace of the precedent splitting (useful in BSPT.cpp until PPT)
	doubleVector split_offset;

	LinearSystem();
	LinearSystem(vector< vector<double> > A, vector< double > b);
	LinearSystem(lst vars, lst constraints);

	vector< vector<double> > getA(); //return A
	vector<double> getb(); 			//return b
	double getA(int i, int j); 		//return A[i][j]
	double getb(int i); 			//return b[i]
	void setb(vector<double> offset){this->b = offset;}

	double minLinearSystem(lst vars, ex obj_fun);
	double maxLinearSystem(lst vars, ex obj_fun);
	double maxLinearSystem(vector< double > obj_fun_coeffs);
	bool isEmpty();
	bool isEmpty_temp(); // TEMP (is_Empty is not working as expected)
	LinearSystem* appendLinearSystem(LinearSystem *LS);
	vector<bool> redundantCons();

	int dim(){ if(!this->isEmpty()){ return this->A[0].size(); }else{ return 0;}	};
	int size(){ return this->b.size(); };

	double volBoundingBox();

	void print();
	std::string print_string();
	void plotRegion();
	void plotRegionToFile(char *file_name, char color);
	void plotRegionT(double t);
	void plotRegion(vector<int> rows, vector<int> cols);

	virtual ~LinearSystem();
};

#endif /* LINEARSYSTEM_H_ */
