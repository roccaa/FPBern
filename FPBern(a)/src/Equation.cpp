/*
 * Equation.cpp
 *
 *  Created on: 15 avr. 2015
 *      Author: alex
 */

#include "Equation.h"


Equation::Equation()
{}

Equation::Equation(	int dim_global, lst vars, lst params, ex eq,
					bool exact, lst q_s, lst beta_s){
	this->base_vars = q_s;
	this->ampl_vars = beta_s;

	this->coeff_bern_numer = NULL;
	this->coeff_bern_denom = NULL;

	this->coeff_bern_max = NULL;
	this->coeff_bern_min = NULL;
	this->exact = exact;
	this->use_guess = false;
	this->dim_global = dim_global;
	this->equationEx = eq;
	this->vars = vars;
	this->params = params;
	lst tempVars;
	intVector varTable;
	this->type = Polynomial_ONLY;
		
	int tempDim = 0;
	for(int k=0;k<vars.nops();k++)
	{
		if(eq.has(vars[k]))
		{
			tempVars.append(vars[k]);
			tempDim++;
			varTable.push_back(k);
		}
	}
	this->varAssociationTable.push_back(varTable);
	//cout << "In equation Constructor\n";
	//cout << "the variables are " << tempVars << endl;
	BaseConverter bc = BaseConverter(tempVars,eq);
	if(tempVars.nops()>1)
	{
		if(this->exact == false){
			this->explicitBC_coeff = bc.getBernCoeffsMatrix();
		}	
		else{
			this->explicitBC_coeff = bc.getBernCoeffs();
		}	
	}
	else
	{
		this->explicitBC_coeff = bc.getBernCoeffs();
	}
	assert(this->varAssociationTable.size() == 1);
	//cout << "Number of bernstein coeff: " << this->explicitBC_coeff.nops() << endl; 	
}


///////////////////////////////////////////////////////////////////////////////////////////////////
Equation::Equation(	int dim_global, lst vars, lst params, ex numer, ex denom,
					bool exact, lst q_s, lst beta_s){
	this->base_vars = q_s;
	this->ampl_vars = beta_s;

	this->coeff_bern_numer = NULL;
	this->coeff_bern_denom = NULL;

	this->coeff_bern_max = NULL;
	this->coeff_bern_min = NULL;

	this->use_guess = false;
	this->dim_global = dim_global;
	this->equationEx = numer/denom;
	this->vars = vars;
	this->params = params;
	this->exact = exact;
	//ex temp = eq.normal();
	lst tempVars;
	int tempDim = 0;
	intVector varTable;
	for(int k=0;k<vars.nops();k++){
		if(this->equationEx.has(vars[k])){
			tempVars.append(vars[k]);
			tempDim++;
			varTable.push_back(k);
		}
	}
	if(varTable.size()!=0) this->varAssociationTable.push_back(varTable);
	//cout << "computing numer denom\n"; 
	//ex nd = eq.numer_denom();
	BaseConverter bc = BaseConverter(tempVars,numer, denom);
	this->explicitBC_coeff = bc.getRationalBernCoeffs();
	//cout << "Nb bernCoeffs = " << this->explicitBC_coeff.nops() << endl; 
	this->type = RationalPolynomial_ONLY;

}

pair<double,double> Equation::optimize(Parallelotope *set, LinearSystemSet *paramSet)
{
	//cout << "In [optimize()]" << endl;
	bern_info optimization;
	optimization = optimize_explicit(set, paramSet);
	//cout << "recup results\n";
	pair<double,double> res;
	if(this->exact==true){
		res.first = ex_to<numeric>(optimization.min_e).to_double();
		res.second = ex_to<numeric>(optimization.max_e).to_double();
	}
	else{
		res.first = optimization.min;
		res.second = optimization.max;
	}
	
	//cout << "delete coeffs list\n";
	delete(optimization.coeff_max);
	delete(optimization.coeff_min);

	return res;
}


pair<ex,ex> Equation::optimize_exact(Parallelotope *set, LinearSystemSet *paramSet)
{
	//cout << "In [optimize_exact()]" << endl;
	bern_info optimization;

	//cout << "entering optimization code\n";
	optimization = optimize_explicit(set, paramSet);

	//cout << "recup results\n";
	pair<ex,ex> res;
	res.first = optimization.min_e;
	res.second = optimization.max_e;
	//cout << "delete coeffs list\n";
	//delete(optimization.coeff_max);
	//delete(optimization.coeff_min);

	return res;
}

double abs_sum(ex p)
{
	//cout << "in abs_sum func \n"; 
	ex poly = p.expand();
	int nbmonome = poly.nops();
	double sum = 0;
	//ex exsum = 0;
	for(int i=0;i<nbmonome;i++)
	{
		ex m = poly.op(i);
		if(m.nops() < 2)
		{
			if( is_a<numeric>( (m) ) )
			{
				sum = sum + abs(ex_to<numeric>((m)).to_double());
				//exsum = exsum + abs(poly.op(i));
			}
			else
			{
				sum = sum + 1;
				//exsum = exsum + 1;
			}
		}
		else
		{
			sum = sum + abs(ex_to<numeric>((m).op(1)).to_double());
			//exsum = exsum + abs((poly.op(i)).op(1));
		}

	}

    //cout << "exsum = " << ex_to<numeric>(exsum).to_double() << endl;
	//cout << "sum = " << sum << endl;
	//assert(0);
	return sum;
	//return ex_to<numeric>(exsum).to_double();
}


ex exact_abs_sum(ex p)
{
	//cout << "in exact_abs_sum func \n"; 
	ex poly = p.expand();
	int nbmonome = poly.nops();
	ex sum = 0;
	//ex exsum = 0;
	for(int i=0;i<nbmonome;i++)
	{
		ex m = poly.op(i);
		if(m.nops() < 2)
		{
			if( is_a<numeric>( (m) ) )
			{
				sum = sum + abs(ex_to<numeric>((m)));
				//exsum = exsum + abs(poly.op(i));
			}
			else
			{
				sum = sum + 1;
				//exsum = exsum + 1;
			}
		}
		else
		{
			sum = sum + abs(ex_to<numeric>((m).op(1)));
			//exsum = exsum + abs((poly.op(i)).op(1));
		}

	}

    //cout << "exsum = " << ex_to<numeric>(exsum).to_double() << endl;
	//cout << "sum = " << sum << endl;
	//assert(0);
	return sum;
	//return ex_to<numeric>(exsum).to_double();
}


bern_info Equation::optimize_explicit(Parallelotope *set, LinearSystemSet *paramSet)
{
	bern_info result;

	lst coeffs_symb = this->explicitBC_coeff;
	//cout << "BC symb = \n " << coeffs_symb << "\n";
	result.partial = false;
	result.method = method;

	result.coeff_max = new doubleVector(coeffs_symb.nops());
	result.coeff_min = new doubleVector(coeffs_symb.nops());

	lst subs_map;
	poly_values p;

	int type;
	if(set != NULL) // If the set is not the unit box (not called from optimize(paramSet) function
	{
		p = set->get_poly_values();
		for(int i=0;i<p.lenghts.size();i++){
			subs_map.append(this->base_vars[i] == p.base_vertex[i]);
			subs_map.append(this->ampl_vars[i] == p.lenghts[i]);
		}
		//cout << "subs_map = \n" << subs_map << endl;
		type = set->type;
	}
	type = PARALLELOTOPE;
	//else
	//	type = BOX;


	double maximum;
	double minimum;
	
	switch(type)
	{
		case BOX: // Box only the UnitBox --> Non-Unit Box not Handled
			maximum = -INFINITY;
			minimum =  INFINITY;
			for(int j=0; j<coeffs_symb.nops();j++)
			{
				coeffs_symb[j] = coeffs_symb[j].subs(subs_map,subs_options::no_pattern);
				result.coeff_max->at(j) = optimize_max(coeffs_symb[j],paramSet,this->params);
				result.coeff_min->at(j) = optimize_min(coeffs_symb[j],paramSet,this->params);

				minimum  = min(result.coeff_min->at(j),minimum);
				maximum  = max(result.coeff_max->at(j),maximum);
			}

			result.min = minimum;
			result.max = maximum;
			break;
		default:
			maximum = -INFINITY;
			minimum = INFINITY; // -INFINITY
			//cout << "setting initial values for max and min\n"; 
			parser p;
			ex minimum_e; 
			if(exact ==true )
				minimum_e = p("2^128"); // Max value in Ginac ?
			ex maximum_e = -minimum_e;
			int nbcoeff = coeffs_symb.nops();
			//cout << " In optimize function\n";
			//cout << "Number of bernstein coeff = " << nbcoeff << endl;
			//for(int j=0; j<coeffs_symb.nops();j++)
			double sum_d; 
			ex sum_e; 
			for (lst::const_iterator c = coeffs_symb.begin(); c != coeffs_symb.end(); ++c)
			//for(int j=0; j<nbcoeff;j++)
			{
				//cout << "substitution\n";
				//coeffs_symb[j] = coeffs_symb[j].subs(subs_map,subs_options::no_pattern);
				//cout << "Linear optimization\n"; 
				//cout << "#######################################" << endl;
				//cout << "the current coeff is:\n";
				//cout << (*c) << endl;
				//cout << "The parameters are: \n";
				//cout << this->params << endl;
				//result.coeff_max->at(j) = optimize_max(coeffs_symb[j],paramSet,this->params);
				//result.coeff_min->at(j) = optimize_max(-coeffs_symb[j],paramSet,this->params);
				//minimum  = max(result.coeff_min->at(j),minimum);
				//maximum  = max(result.coeff_max->at(j),maximum);
					
				if(this->exact==true){
					sum_e = exact_abs_sum(*c);
					minimum_e  = min(-sum_e,minimum_e);
					maximum_e  = max(sum_e,maximum_e);
				}
				else{
					sum_d = abs_sum(*c);
					minimum  = min(-sum_d,minimum);
					maximum  = max(sum_d,maximum);
					//cout  <<"minimum = " << minimum <<"\n";
					//cout << "maximum = " << maximum << "\n";
					//cout << "#######################################" << endl;
				}
				//sum = abs_sum(coeffs_symb[j]);

				
				
			}
			//cout << "done\n";
			result.min = minimum;
			result.max = maximum;
			result.min_e = minimum_e;
			result.max_e = maximum_e;
			break;
	}

	return result;
}

double Equation::max_coeff(lst coeff_list, Polyhedron *set, LinearSystemSet * paramSet)
{
	double maximum = -INFINITY;
	lst coeff_s = coeff_list;
	doubleVector coeff = doubleVector(coeff_list.nops());
	lst subs_map;
	poly_values p = set->get_poly_values();
	vector<double> q = vector<double>(p.base_vertex.size());
	vector<double> beta = vector<double>(p.lenghts.size());
	for(int i=0;i<p.lenghts.size();i++){
		subs_map.append(this->base_vars[i] == p.base_vertex[i]);
		subs_map.append(this->ampl_vars[i] == p.lenghts[i]);
	}
	for(int j=0; j<coeff_s.nops();j++) // TODO
	{
		coeff_s[j] = coeff_s[j].subs(subs_map,subs_options::no_pattern);
		coeff[j] = optimize_max(coeff_s[j],paramSet,this->params);
		maximum  = max(coeff[j],maximum);

	}
	return maximum;
}

double Equation::min_coeff(lst coeff_list, Polyhedron *set, LinearSystemSet * paramSet)
{
	double minimum = INFINITY;
	lst coeff_s = coeff_list;
	doubleVector coeff = doubleVector(coeff_list.nops());
	lst subs_map;
	poly_values p = set->get_poly_values();
	vector<double> q = vector<double>(p.base_vertex.size());
	vector<double> beta = vector<double>(p.lenghts.size());
	for(int i=0;i<p.lenghts.size();i++){
		subs_map.append(this->base_vars[i] == p.base_vertex[i]);
		subs_map.append(this->ampl_vars[i] == p.lenghts[i]);
	}
	for(int j=0; j<coeff_s.nops();j++) // TODO
	{
		coeff_s[j] = coeff_s[j].subs(subs_map,subs_options::no_pattern);
		coeff[j] = optimize_min(coeff_s[j],paramSet,this->params);
		minimum  = min(coeff[j],minimum);

	}
	return minimum;
}

double Equation::optimize_max(ex num_cs,LinearSystemSet *paramSet,lst params)
{
	if(is_a<numeric>(num_cs))
	{
		return ex_to<numeric>(num_cs).to_double();
	}
	else
	{
		return paramSet->getSet()[0]->maxLinearSystem(this->params,num_cs);
	}
}

double Equation::optimize_min(ex num_cs,LinearSystemSet *paramSet,lst params)
{
	if(is_a<numeric>(num_cs))
	{
		return ex_to<numeric>(num_cs).to_double();
	}
	else
	{
		return paramSet->getSet()[0]->minLinearSystem(this->params,num_cs);
	}
}


void Equation::print_ex(void)
{
	cout << this->equationEx << endl;
}




Equation::~Equation()
{


}
