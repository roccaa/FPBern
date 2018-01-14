/*
 * Equation.h
 *
 *  Created on: 15 avr. 2015
 *      Author: alex
 */

#ifndef EQUATION_H_
#define EQUATION_H_

#include "LinearSystemSet.h"
#include "Parallelotope.h"
#include "BaseConverter.h"

class Equation {

public:
	Equation();
	Equation(int dim_global,lst vars, lst params, ex eq, bool rationnal, bernstein_computation_method method_BC, lst q_s, lst beta_s);
	~Equation();

	//TODO
	// compute numerical bernstein coeff(poly_value,double interval)
	// compute berstein coeff of a bisection
	// recompute num_berncoeff => after split length and base vertex change

	pair<double,double> optimize(Parallelotope *set, LinearSystemSet *paramSet);

	bern_info optimize_explicit(Parallelotope *set, LinearSystemSet *paramSet);

	double optimize_max(ex num_cs,LinearSystemSet *paramSet,lst params);
	double optimize_min(ex num_cs,LinearSystemSet *paramSet,lst params);
	double max_coeff(lst coeff_list, Polyhedron *set, LinearSystemSet * paramSet);
	double min_coeff(lst coeff_list, Polyhedron *set, LinearSystemSet * paramSet);

	void constructor_explicit_sym(int dim_global,lst vars, lst params, ex eq, bool fully_rationnal);

	void print_ex(void);

	double getMaxBound(void){return this->maxBound;};
	double getMinBound(void){return this->minBound;};
	bool use_guess;

	vector<ex>* coeff_bern_numer;
	vector<ex>* coeff_bern_denom;

	vector<ex>* coeff_bern_max;
	vector<ex>* coeff_bern_min;


private:
	lst explicitBC_num;
	lst explicitBC_denom;

	double maxBound;
	double minBound;
	int dim_global;
	lst vars;
	lst params;
	ex equationEx; //Ginac expression
	vector<intVector> varAssociationTable;
	equation_type type;
	bernstein_computation_method method;
	lst base_vars;
	lst ampl_vars;
	lst explicitBC_coeff;


};




#endif /* EQUATION_H_ */
