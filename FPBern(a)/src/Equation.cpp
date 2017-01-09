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
					bernstein_computation_method method_BC,
					lst q_s, lst beta_s){
	this->method = method_BC;
	this->base_vars = q_s;
	this->ampl_vars = beta_s;

	this->coeff_bern_numer = NULL;
	this->coeff_bern_denom = NULL;

	this->coeff_bern_max = NULL;
	this->coeff_bern_min = NULL;

	switch(method_BC)
	{
		case IMPLICIT_NUM:
			cout << "NO method for this tool"<<endl;
			break;
		case IMPLICIT_SYM:
			cout << "NO method for this tool"<<endl;
			break;
		case EXPLICIT_SYM:
			constructor_explicit_sym(dim_global, vars, params, eq, false);
			break;
		default:
			constructor_explicit_sym(dim_global, vars, params, eq, false);
			break;
	}
}

void Equation::constructor_explicit_sym(int dim_global,lst vars, lst params, ex eq, bool fully_rationnal)
{
	this->use_guess = false;
	this->dim_global = dim_global;
	this->equationEx = eq;
	this->vars = vars;
	this->params = params;
	ex tempEq = eq.expand();
	lst tempVars;
	intVector varTable;
	this->type = Polynomial_ONLY;
	if(!fully_rationnal){ // Just polynomial
		int tempDim = 0;
		for(int k=0;k<vars.nops();k++)
		{
			if(tempEq.has(vars[k]))
			{
				tempVars.append(vars[k]);
				tempDim++;
				varTable.push_back(k);
			}
		}
		this->varAssociationTable.push_back(varTable);
		//cout << "In equation Constructor\n";
		//cout << "the variables are " << tempVars << endl;
		BaseConverter bc = BaseConverter(tempVars,tempEq);
		if(tempVars.nops()>1)
		{
			this->explicitBC_coeff = bc.getBernCoeffsMatrix();
		}
		else
		{
			this->explicitBC_coeff = bc.getBernCoeffs();
		}
		assert(this->varAssociationTable.size() == 1);
		//cout << "Number of bernstein coeff: " << this->explicitBC_coeff.nops() << endl; 		
	}
	else{
		cout << " \n ####### WARNING !! ####\n /!\ Execution terminating !\n";
		cout << "====>>> Rationnal Function Cannot be computed with this method \n";
		cout << "====>>> getBernCoeffsMatrix, and getBernCoeffs may not give the same number of coeff for numerator and denominator\n";
		cout << "====>>> getBernCoeffsMatrix, and getBernCoeffs need a function to compute higher degrees bernstein decomposition\n";
		assert(0);
		this->type = RationalPolynomial_ONLY;
	}
}

pair<double,double> Equation::optimize(Parallelotope *set, LinearSystemSet *paramSet)
{

	bern_info optimization;
	switch(method)
	{
		case IMPLICIT_NUM:
			break;

		case IMPLICIT_SYM:
			break;

		case EXPLICIT_SYM:
			optimization = optimize_explicit(set, paramSet);
			break;

		default:
			break;
	}
	pair<double,double> res;
	res.first = optimization.min;
	res.second = optimization.max;

	delete(optimization.coeff_max);
	delete(optimization.coeff_min);

	if(this->use_guess && optimization.min_coeffs_guess != NULL)
	{
		delete(optimization.min_coeffs_guess);
	}
	if(this->type != RationalPolynomial_ONLY && method == IMPLICIT_SYM && optimization.max_coeffs_guess != NULL)
		delete(optimization.max_coeffs_guess);
	if(this->type == RationalPolynomial_ONLY && optimization.coeff_numer != NULL)
	{
		delete(optimization.coeff_numer);
		delete(optimization.coeff_denom);
	}

	return res;
}

double abs_sum(ex p)
{
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


bern_info Equation::optimize_explicit(Parallelotope *set, LinearSystemSet *paramSet)
{
	bern_info result;

	lst coeffs_symb = this->explicitBC_coeff;
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
		cout << "subs_map = \n" << subs_map << endl;
		type = set->type;
	}
	type = PARALLELOTOPE;
	//else
	//	type = BOX;


	double maximum;
	double minimum;
	// Remainder-----> Do not Handle rational function with this method !
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
			int nbcoeff = coeffs_symb.nops();
			//cout << " In optimize function\n";
			//cout << "Number of bernstein coeff = " << nbcoeff << endl;
			//for(int j=0; j<coeffs_symb.nops();j++)
			double sum;  
			for (lst::const_iterator c = coeffs_symb.begin(); c != coeffs_symb.end(); ++c)
			//for(int j=0; j<nbcoeff;j++)
			{
				//cout << "substitution\n";
				//coeffs_symb[j] = coeffs_symb[j].subs(subs_map,subs_options::no_pattern);
				//cout << "Linear optimization\n"; 
				//cout << "the current coeff is:\n";
				//cout << coeffs_symb[j] << endl;
				//cout << "The parameters are: \n";
				//cout << this->params << endl;
				//result.coeff_max->at(j) = optimize_max(coeffs_symb[j],paramSet,this->params);
				//result.coeff_min->at(j) = optimize_max(-coeffs_symb[j],paramSet,this->params);
				//minimum  = max(result.coeff_min->at(j),minimum);
				//maximum  = max(result.coeff_max->at(j),maximum);
				
				sum = abs_sum(*c);
				//sum = abs_sum(coeffs_symb[j]);
				minimum  = min(-sum,minimum);
				maximum  = max(sum,maximum);
				
				
			}
			//cout << "done\n";
			result.min = minimum;
			result.max = maximum;
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
