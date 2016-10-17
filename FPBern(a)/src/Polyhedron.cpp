/*
 * Polyhedron.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: dreossi
 */

#include "Polyhedron.h"

Polyhedron::Polyhedron() {
	// TODO Auto-generated constructor stub
	this->type = 2;
	this->re_used = false;
}

// Get functions
lst Polyhedron::getQ(){ return this->vars[0]; }
lst Polyhedron::getAlpha(){ return this->vars[1]; }
lst Polyhedron::getBeta(){  return this->vars[2]; }
lst Polyhedron::getGeneratorFunction(){	return this->generator_function; }
int Polyhedron::getDim(){ return this->dim; }
vector< vector< double > > Polyhedron::getTemplate(){ return this->template_matrix; }
void Polyhedron::add_BernsteinInfo(BernCoeff_for_Box * info){info_BernsteinCoeff.push_back(info);};
vector<BernCoeff_for_Box *> *  Polyhedron::get_bernsteinInfo(){return &info_BernsteinCoeff;}
/*Box * Polyhedron::compute_reachedBox(void)
{
	assert(this->info_BernsteinCoeff.empty()==false);
	doubleVector bl = doubleVector(this->dim);
	doubleVector bu = doubleVector(this->dim);
	assert(this->info_BernsteinCoeff.size() == dim); // Meme nombre dynamiques et de variables ...
	for(int i=0; i<this->info_BernsteinCoeff.size(); i++){
		bl[i] = this->info_BernsteinCoeff[i]->min;
		bu[i] = this->info_BernsteinCoeff[i]->max;
	}
	return new Box(bl,bu);
}*/
void Polyhedron::clear_BCinfo(void)
{
	//cout << "clearBCinfo\n";
	//cout << "size BCinfo = " << info_BernsteinCoeff.size() << endl;
	for(int i=0;i<info_BernsteinCoeff.size();i++){

		if(info_BernsteinCoeff[i]->num_coeffmin != NULL)
		{
			info_BernsteinCoeff[i]->num_coeffmin->clear();
			delete(info_BernsteinCoeff[i]->num_coeffmin);

		}
		if(info_BernsteinCoeff[i]->num_coeffmax != NULL)
		{
			info_BernsteinCoeff[i]->num_coeffmax->clear();
			delete(info_BernsteinCoeff[i]->num_coeffmax);

		}
		if(info_BernsteinCoeff[i]->denom_coeffmin != NULL)
		{
			info_BernsteinCoeff[i]->denom_coeffmin->clear();
			delete(info_BernsteinCoeff[i]->denom_coeffmin);
		}
		if(info_BernsteinCoeff[i]->denom_coeffmax != NULL)
		{
			info_BernsteinCoeff[i]->denom_coeffmax->clear();
			delete(info_BernsteinCoeff[i]->denom_coeffmax);
		}
	}
	//cout <<"inside info deleted\n";
	while(!this->info_BernsteinCoeff.empty()){
		delete(this->info_BernsteinCoeff[this->info_BernsteinCoeff.size()-1]);
		this->info_BernsteinCoeff.pop_back();
	}
	//cout <<"container deleted\n";
	info_BernsteinCoeff.clear();
	//cout <<"done\n";
}


Polyhedron::~Polyhedron() {
	// TODO Auto-generated destructor stub
}

