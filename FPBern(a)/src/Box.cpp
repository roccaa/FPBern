/*
 * Box.cpp
 *
 *  Created on: Nov 3, 2014
 *      Author: dreossi
 */

#include "Box.h"


Box::Box(vector<double> l,vector<double> u){
	this->l = l;
	this->u = u;
	bounds2const(l,u);
	this->dim = l.size();
	this->motherLocId = -2;
	this->LocId = -2;
	this->re_used = false;
	this->type = BOX;
	//cout << "new_box !\n";
}

Box::Box(LinearSystem* linSystem){
	this->linSystem = linSystem;
	cout << "generating the bounds \n";
	linSystem->print();
	vector< vector<double> > tempBounds = this->const2bounds(this->linSystem);
	cout << "Box generated \n";
	this->l = tempBounds[0];
	this->u = tempBounds[1];
	this->dim = l.size();
	this->motherLocId = -2;
	this->LocId = -2;
	this->re_used = false;
	this->type = BOX;

}

Box::Box(vector<lst> vars) {

	if(vars.size() != 3){
		cout<<"Box::Box : vars must contain 3 collections of variable names (q,alpha,beta)";
		exit (EXIT_FAILURE);
	}

	this->vars.push_back(vars[0]);
	this->vars.push_back(vars[1]);
	this->vars.push_back(vars[2]);

	// get the dimension of the box
	this->dim = vars[0].nops();

	// and store its variable names
	for(int i=0; i<3; i++){
		if((signed)vars[i].nops() != this->dim){
			cout<<"Box::Box : vars["<<i<<"] must have "<<this->dim<<" variables";
			exit (EXIT_FAILURE);
		}
	}

	// extract variable names
	lst q = this->vars[0];
	lst alpha = this->vars[1];
	lst beta = this->vars[2];

	// initialize generator function
	for(int i=0; i<this->dim; i++){
		this->generator_function.append(q[i]);
	}

	// create the generation function acclumulatin the versor values
	for(int i=0; i<this->dim; i++){
		this->generator_function[i] = this->generator_function[i] + alpha[i]*beta[i];
	}
	this->re_used = false;
	this->type = BOX;
}

// convert from generator to constraints representation
// q : numeric base vertex
LinearSystem* Box::gen2const(vector<double> q, vector<double> beta){

	if((signed)q.size() != this->dim){
		cout<<"Box::gen2const : q must have dimension "<<this->dim;
		exit (EXIT_FAILURE);
	}

	vector< vector<double> > Lambda;
	vector< double > d;

	// initialize template Lambda
	for(int i=0; i<this->dim; i++){
		vector<double> Lambda_i (this->dim,0);
		Lambda_i[i] = 1;
		Lambda.push_back(Lambda_i);
		d.push_back(q[i]+beta[i]);
	}
	for(int i=0; i<this->dim; i++){
		vector<double> Lambda_i (this->dim,0);
		Lambda_i[i] = -1;
		Lambda.push_back(Lambda_i);
		d.push_back(-q[i]);
	}

	LinearSystem *LS = new LinearSystem(Lambda,d);
	return LS;
}

poly_values Box::const2gen(LinearSystem *constr){

	poly_values res;
	return res;

}


vector< vector<double> > Box::const2bounds(LinearSystem* linSystem){
	cout << "BOX::const2bounds\n";
	vector< vector<double> > result;
	int n = linSystem->getb().size();
	vector<double> b = linSystem->getb();
	cout << "b = " << b;
	vector<double> bl = vector<double>(n/2);
	vector<double> bu = vector<double>(n/2);
	cout << "n/2 = " << n/2 << endl;
	for(int i=0;i<n;i++){
		int c = i%(n/2);
		if(i<n/2)
			bl[i] = -b[i];
		if(i>=n/2)
			bu[c] = b[i];
		cout << "b["<<i<<"] = " << b[i] << endl;
	}
	result.push_back(bl);result.push_back(bu);

	return result;
}

LinearSystem* Box::bounds2const(vector<double> l, vector<double> u){

	int n = l.size();
	vector< vector<double> > A;
	vector<double> b = vector<double>(2*n);
	for(int i=0;i<n;i++){
		vector<double> Ai = vector<double>(n,0.0);
		Ai[i] = -1;
		A.push_back(Ai);
	}
	for(int i=0;i<n;i++){
		vector<double> Ai = vector<double>(n,0.0);
		Ai[i] = 1;
		A.push_back(Ai);
	}
	for(int i=0;i<2*n;i++){
		if(i<n)
			b[i] = -l[i];
		else
			b[i] = u[i%n];
	}

	LinearSystem *result = new LinearSystem(A,b);
	this->linSystem = result;

	return result;

}

double Box::get_fullBoxSize(){

	double result = -INFINITY;
	for(int i=0;i<this->l.size();i++){
		if((u[i]-l[i])>result)
			result = u[i]-l[i];
	}
	return result;
}

int Box::compute_MaxSizeDirection(){
	int dir = 0;
	double result = -INFINITY;
	for(int i=0;i<this->l.size();i++){
		if((u[i]-l[i])>result){
			result = u[i]-l[i];
			dir = i;
		}
	}
	return dir;
}

Box* Box::compute_reachedBox(void){
	assert(this->info_BernsteinCoeff.empty()==false);
	doubleVector bl = doubleVector(this->info_BernsteinCoeff.size());
	doubleVector bu = doubleVector(this->info_BernsteinCoeff.size());
	//assert(this->info_BernsteinCoeff.size() == dim); // Meme nombre dynamiques et de variables ...
	for(int i=0; i<this->info_BernsteinCoeff.size(); i++){
		bl[i] = this->info_BernsteinCoeff[i]->min;
		bu[i] = this->info_BernsteinCoeff[i]->max;
	}
//	cout << "reached before fusion Bl= " << bl;
//	cout << "reached before fusion Bu= " << bu;
	return new Box(bl,bu);
}/*
void Box::clear_BCinfo(void){
	for(int i=0;i<info_BernsteinCoeff.size();i++){
		info_BernsteinCoeff[i]->num_coeffmin->clear();
		delete(info_BernsteinCoeff[i]->num_coeffmin);
		info_BernsteinCoeff[i]->num_coeffmax->clear();
		delete(info_BernsteinCoeff[i]->num_coeffmax);
		info_BernsteinCoeff[i]->denom_coeffmin->clear();
		delete(info_BernsteinCoeff[i]->denom_coeffmin);
		info_BernsteinCoeff[i]->denom_coeffmax->clear();
		delete(info_BernsteinCoeff[i]->denom_coeffmax);
	}
	while(!this->info_BernsteinCoeff.empty()){
		delete(this->info_BernsteinCoeff[this->info_BernsteinCoeff.size()-1]);
		this->info_BernsteinCoeff.pop_back();
	}
	info_BernsteinCoeff.clear();
}*/

Box::~Box() {
	// TODO Auto-generated destructor stub
	//cout << "delete BC info in Box\n";
	clear_BCinfo();
	//cout << "delete linSystem in Box\n";
	delete(this->linSystem);
	//cout << "deleting done in Box\n";
}


