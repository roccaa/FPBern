/*
 * Polyhedron.h
 *
 *  Created on: Oct 30, 2014
 *      Author: dreossi
 */

#ifndef POLYHEDRON_H_
#define POLYHEDRON_H_

#include "Common.h"
#include "LinearSystem.h"

struct BernCoeff_for_Box{
	equation_type type;
	vector<double> * num_coeffmin;
	vector<double> * num_coeffmax;
	vector<double> * denom_coeffmin;
	vector<double> * denom_coeffmax;
	double min;
	double max;

};


class Polyhedron {

protected:

	int dim;							// dimension of the parallelotope
	vector<lst> vars;					// variables appearing in generato function
										// vars[0] q: base vertex
										// vars[1] alpha : free variables \in [0,1]
										// vars[2] beta : generator amplitudes
	lst generator_function;				// generator function
	vector< vector<double> > u;			// versors
	vector< vector<double> > template_matrix;	// Template matrix
	LinearSystem * linSystem;
	poly_values num_genDescription;
	int index_associatedEqSyst;
	lst directions_functions;

public:

	Polyhedron();

	bool re_used;
	intVector previsous_set_loc;

	lst getGeneratorFunction();
	lst getQ();
	lst getAlpha();
	lst getBeta();

	vector<lst> getVars(){return this->vars;}
	lst get_dir_f(){return this->directions_functions;}
	int getAssociatedEqSyst(){return this->index_associatedEqSyst;}
	void setAssociatedEqSyst(int index){this->index_associatedEqSyst = index;}
	vector<doubleVector> getU(){return this->u;}
	void reset_linsys(){if(this->linSystem!=NULL){delete(this->linSystem);this->linSystem=NULL;}}
	int type;
	int getDim();
	vector< vector< double > > getTemplate();
	LinearSystem * getLinSyst(){return linSystem;};
	virtual LinearSystem* gen2const(vector<double> q, vector<double> beta){ return 0; }
	virtual poly_values const2gen(LinearSystem *constr){ poly_values p; return p; };
//	double compute_directionSize(int direction){assert(0);};
	virtual int getMaxSize(){ return 0; };
	virtual vector<ex> getConvCombs(int i){ vector<ex> c; return c; };

	vector<BernCoeff_for_Box *> info_BernsteinCoeff;
	void add_BernsteinInfo(BernCoeff_for_Box * info);
	vector<BernCoeff_for_Box *> *  get_bernsteinInfo();
	void clear_BCinfo(void);
	Polyhedron *compute_reachedBox(void){cout << "should Not go there\n";};
	poly_values get_poly_values(){return this->num_genDescription;}


	virtual ~Polyhedron();
};

#endif /* POLYHEDRON_H_ */
