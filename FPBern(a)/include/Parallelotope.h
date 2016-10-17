/*
 * Parallelotope.h
 *
 *  Created on: Oct 30, 2014
 *      Author: dreossi
 */

#ifndef PARALLELOTOPE_H_
#define PARALLELOTOPE_H_

#include "Box.h"
#include "LinearSystem.h"

class Parallelotope : public Polyhedron {

private:

	vector<double> hyperplaneThroughPts(vector< vector<double> > pts);	// find hyper plane passing through pts
	vector<double> lst2vec(ex list);									// convert a lst to a vector of doubles
	double euclidNorm(vector<double> v);

	vector< double > actual_base_vertex;
	vector< double > actual_lenghts;							// compute euclidean norm

	vector<lst> vars;					// variables appearing in generato function
											// vars[0] q: base vertex
											// vars[1] alpha : free variables \in [0,1]
											// vars[2] beta : generator amplitudes
	lst generator_function;				// generator function
	vector< vector<double> > u;			// versors
	bool re_used;
	int dim;							// dimension of the parallelotope

	vector< vector<double> > template_matrix;	// Template matrix
	LinearSystem * linSystem;
	poly_values num_genDescription;

public:
	int type;
	lst getQ(){ return this->vars[0]; }
	lst getAlpha(){ return this->vars[1]; }
	lst getBeta(){  return this->vars[2]; }
	lst getGeneratorFunction(){	return this->generator_function; }
	int getDim(){ return this->dim; }
	vector< vector< double > > getTemplate(){ return template_matrix; }
	Parallelotope(vector<lst> vars, vector< vector<double> > u);
	Parallelotope(vector<lst> vars, LinearSystem *LS);

	// Representation conversion
	LinearSystem* gen2const(vector<double> q, vector<double> beta);		// from generator to constraints
	poly_values const2gen(LinearSystem *LS);							// from constraints to generators

	Box * overBox;
	vector_matrix vertices_list;

	//vector<BernCoeff_for_Box *> info_BernsteinCoeff;
	//void add_BernsteinInfo(BernCoeff_for_Box * info){info_BernsteinCoeff.push_back(info);};
	//vector<BernCoeff_for_Box *> *  get_bernsteinInfo(){return &info_BernsteinCoeff;}
	Box *compute_reachedBox(void);
	vector_matrix generate_template(vector_matrix normals);
	void recompute_gen_func(vector_matrix u);
	void set_template(vector_matrix templt){this->template_matrix = templt;}
	void set_linSys(LinearSystem *linSys){delete(this->linSystem);this->linSystem=linSys;}
	void add_polyValues(poly_values p);
	void  add_polyValues(poly_values p, bool do_overBox);
	vector<doubleVector> getU(){return this->u;}
	Box *compute_overBox();
	LinearSystem * getLinSyst(){return linSystem;};
	void clear_BCinfo(void);
	poly_values get_poly_values(){return this->num_genDescription;}

	virtual ~Parallelotope();
};

#endif /* PARALLELOTOPE_H_ */
