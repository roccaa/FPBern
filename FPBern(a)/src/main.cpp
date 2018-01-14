/*
 * test_Implicit.cpp
 * Test l'optimisation polynomial avec la forme implicit de bernstein
 *  Created on: 10 dec 2015
 *      Author: alex
 */


#include <stdio.h>
#include <iostream>
#include <fstream>

#include "Equation.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

using namespace std;


int main(int argc,char** argv){


////////////////////////////////////////////////////////////////////////////// 
/*
	clock_t tStart = clock();	
	possymbol   x1("x1"), x2("x2"), x3("x3"), x4("x4");
	lst vars;
	vars = x1,x2;
	lst params;
	params = x3,x4;

	ex eq = ((3*pow(x1,2)/5)+x2)/((5*x1)/3+pow(x2,2)+2)*x3 + (x1/(2+x2))*x4;
	//ex eq = (3*pow(x1,2)/5+x2)*x3 + (x1)*x4;
	eq.expand();

	vector<double> pAi (params.nops(),0);
	vector< vector<double> > pA (2*params.nops(),pAi);
	vector<double> pb (2*params.nops(),0);
	// left constraint .... <= g(x) --> -g(x) <= -....
	for(int i=0;i<2*params.nops();i++)
	{
		if(i<(params.nops()))
			pA[i][i] = -1;
		else
			pA[i][i-(params.nops())] = 1;

		pb[i] = 1;
	}
	LinearSystem *parameters = new LinearSystem(pA,pb);
	LinearSystemSet *paramSet = new LinearSystemSet(parameters);
	int dim_global = 2	;

	// Calcul with parallelotope description (here a [0;1] box)
	// DEFINITION Of the parallelotopes
	possymbol q1("q1"), q2("q2"), q3("q3"), q4("q4"), q5("q5"), q6("q6"), q7("q7");
	possymbol a1("a1"), a2("a2"), a3("a3"), a4("a4"), a5("q5"), a6("q6"), a7("q7");
	possymbol b1("b1"), b2("b2"), b3("b3"), b4("b4"), b5("q5"), b6("q6"), b7("q7");
	lst qs,as,bs; // Symbolic variables for linear transformation
	qs = q1,q2;
	as = a1,a2;
	bs = b1,b2;
	vector<lst> set_vars;
	set_vars.push_back(qs);
	set_vars.push_back(as);
	set_vars.push_back(bs);
	doubleVector row = doubleVector(2,0);
	vector<doubleVector> u= vector<doubleVector>(2,row);
	u[0][0] = 1;
	u[1][1] = 1;
	Parallelotope* para = new Parallelotope(set_vars,u);
	doubleVector bV = {1,1};
	doubleVector ampl = {2,2};

	poly_values p;
	p.base_vertex = bV;
	p.versors = u;
	p.lenghts = ampl;
	para->add_polyValues(p);
	
	cout << "Init set built\n";

	// Computing f°g
	lst generator_function = para->getGeneratorFunction();
	cout << "Getting Generator functions\n";
	lst sub;
	ex fog;
	for(int i=0; i<(signed)generator_function.nops(); i++){
		sub.append(vars[i] == generator_function[i]);
	}
	fog = eq.subs(sub,subs_options::no_pattern);
	lst subs_map;
	for(int i=0;i<p.lenghts.size();i++){
		subs_map.append(qs[i] == p.base_vertex[i]);
		subs_map.append(bs[i] == p.lenghts[i]);
	}
	//cout << "subs_map = \n" << subs_map << endl;
	fog =  fog.subs(subs_map,subs_options::no_pattern);
	fog =  fog.expand();

	cout << "Composition with generator\n";
	bool rational = true;
	Equation* eq_p = new Equation(dim_global, as, params, fog, rational, EXPLICIT_SYM, qs, bs);
	cout << "Equation created:\n";
	cout << fog << " \n" ;
	tStart = clock();
	cout << "Opti\n";
	pair<double,double>  res_explicit = eq_p->optimize(para,paramSet);
	cout << "Time method explicit Parallelotopes: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
	cout << "Result from optimization of: \n" << eq << endl;
	cout << "min = " << -1*res_explicit.first << endl;
	cout << "max = " << res_explicit.second << endl;
	
	cout << "segfault ?" << endl;
    //assert(0);
*/
 /////////////////////////////////////////////////////////////////////////////////  
	if(argc == 1)
	{
		cout << " No input files -- please give a .ini input file with a polynomial program to start\n";
		cout << ">>>>>> shutting down ...";
		return -1;
	}

	for(int i=1;i<argc;i++)
	{
		std::string file = argv[i];
		boost::property_tree::ptree pt;
	    boost::property_tree::ini_parser::read_ini((file+".ini"), pt);
	   // cout << "Begin reading\n";
	    string name =  pt.get<string>("OPTIONS.name");
	    int exp = pt.get<int>("OPTIONS.precision");
	    double err = std::pow(2,-exp);
	    //cout << "err = " << err << endl;
	    int dim_vars = pt.get<int>("OPTIONS.nbvars");
	    int dim_param = pt.get<int>("OPTIONS.nberrors");
	    symtab table_symbols;

	    lst vars;
	    lst params;
		for(int j=1;j<=dim_vars+dim_param;j++)
		{
			string s = "x"+boost::lexical_cast<std::string>(j);
			symbol x(s);
			if(j<=dim_vars)
				vars.append(x);
			else
				params.append(x);
		    table_symbols[s] = x;
		}

	    string p_string  = pt.get<string>("PROGRAM.function");
	    //cout << p_string << endl;
	    parser reader(table_symbols);
	    ex p = reader(p_string);
		bool rational = false;
		if(p.denom() != 1){
			rational = true;
		}
		else
			p = p.expand();
		
	    //cout << p << endl;

	    string separateur = ",";
	    string string_bl = pt.get<string>("PROGRAM.input_bl");
	    doubleVector bl = read_doubleVector(&string_bl,separateur,dim_vars);
	    //cout << bl;
	    string string_bu = pt.get<string>("PROGRAM.input_bu");
	    doubleVector bu = read_doubleVector(&string_bu,separateur,dim_vars	);
	    //cout << bu;

	    vector<double> pAj (params.nops(),0);
		vector< vector<double> > pA (2*params.nops(),pAj);
		vector<double> pb (2*params.nops(),0);
		for(int j=0;j<2*params.nops();j++)
		{
			if(j<(params.nops()))
				pA[j][j] = -1;
			else
				pA[j][j-(params.nops())] = 1;

			pb[j] = 1;
		}

		LinearSystem parameters = LinearSystem(pA,pb);
		LinearSystemSet paramSet = LinearSystemSet(&parameters);

		lst qs,as,bs; // Symbolic variables for linear transformation
		for(int j=1;j<=dim_vars;j++)
		{
			string sq = "q"+boost::lexical_cast<std::string>(j);
			string sa = "a"+boost::lexical_cast<std::string>(j);
			string sb = "b"+boost::lexical_cast<std::string>(j);
			symbol q(sq); symbol a(sa); symbol b(sb);
			table_symbols[sq] = q; table_symbols[sa] = a; table_symbols[sb] = b;
			qs.append(q); as.append(a); bs.append(b);
		}
		vector<lst> set_vars;
		set_vars.push_back(qs);
		set_vars.push_back(as);
		set_vars.push_back(bs);
		doubleVector row = doubleVector(dim_vars,0);
		vector<doubleVector> u= vector<doubleVector>(dim_vars,row);
		for(int j=0;j<dim_vars;j++)
			u[j][j] = 1;
		//cout << "Building init set\n";
		Parallelotope para = Parallelotope(set_vars,u);

		doubleVector bV = bl;
	    doubleVector ampl = bu-bl;

	    poly_values pv;
	    pv.base_vertex = bV;
	    pv.versors = u;
	    pv.lenghts = ampl;
	    //cout << "Addind p values \n";
	    para.add_polyValues(pv);
	    //cout << para.get_poly_values().base_vertex;
	    //cout << "done\n";

	    // Computing f°g
		//cout << "Computing f°g\n";
		lst generator_function = para.getGeneratorFunction();

		lst sub;
		ex fog;
		for(int i=0; i<(signed)generator_function.nops(); i++){
			sub.append(vars[i] == generator_function[i]);
		}
		//cout << "Ready for substitution in p\n";
		//cout << "subs = \n" << sub << endl;
		fog = p.subs(sub,subs_options::no_pattern);

		
		lst subs_map;
		for(int i=0;i<pv.lenghts.size();i++){
			subs_map.append(qs[i] == pv.base_vertex[i]);
			subs_map.append(bs[i] == pv.lenghts[i]);
		}
		//cout << "subs_map = \n" << subs_map << endl;
		fog =  fog.subs(subs_map,subs_options::no_pattern);
		//assert(0);
		fog =  fog.expand();
		//cout << fog << "\n";
		//cout << "######\n";
		//cout << fog.numer() << "\n";
		//cout << "######\n";
 		//cout << fog.denom() << "\n";
		
		//cout << "New program:\n" << fog << endl;
		//cout << "Variable change done -- Creating Equation\n";
		clock_t tStart = clock();
		Equation* eq_p = new Equation(	dim_vars, as, params, fog, rational, EXPLICIT_SYM, qs, bs);
		//Equation* eq_p = new Equation(	dim_vars, as, params, p, EXPLICIT_SYM, qs, bs);
		//cout << "Equation Creation Done -- optimizing Equation\n";
		//pair<double,double>  res_explicit = eq_p->optimize(&para,&paramSet);
		pair<double,double>  res_explicit = eq_p->optimize(NULL,&paramSet);
		//cout << "done!\n";
		double time2 = (double)(clock() - tStart)/CLOCKS_PER_SEC;
		cout << "##### " << name << " #####\n";
		cout << "Total Time: " << time2 << endl;
		//cout << "Result from optimization of: \n" << p << endl;
		//cout << "min = " << -1*res_explicit.first*err << endl;
		//cout << "max = " << res_explicit.second*err << endl;
		cout << "roundoff error = " << max(abs(res_explicit.first*err),abs(res_explicit.second*err)) << endl;

	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	//tStart;
//	double err = std::pow(2,-53);
//
//	possymbol   x1("x1"), x2("x2"), x3("x3"), x4("x4"),
//				x5("x5"), x6("x6");
//	symbol  x7("x7"), x8("x8"), x9("x9"), x10("x10"),
//			x11("x11"), x12("x12"), x13("x13"), x14("x14"), x15("x15"),
//			x16("x16"), x17("x17"), x18("x18"), x19("x19"),
//			x20("x20"), x21("x21"),x22("x22"), x23("x23"), x24("x24"), x25("x25"),
//			x26("x26"), x27("x27"),x28("x28"), x29("x29"), x30("x30"), x31("x31"),
//			x32("x32"), x33("x33"), x34("x34"), x35("x35"), x36("x36"), x37("x37"),
//			x38("x38"), x39("x39"), x40("x40"), x41("x41"),x42("x42"), x43("x43"),
//			x44("x44"),x45("x45"), x46("x46"), x47("x47"),x48("x48");
//
//	// Polynome (fixed parameters)
//	ex eq = ((( (-2/1) * x7 + ( - x9 + ( - x13 + ( - x14 + ( - x15 + ( - x16 + ( - x17 + ( - x18 + ( - x19 + ( - x20 + ( - x28 + ( - x36 + ( - x39 + ( - x42 + ( - x45 + ( - x48)))))))))))))))) * x4) * x1 + ((( (2/1) * x7 + (x9 + (x13 + (x15 + (x16 + (x17 + (x18 + (x19 + (x20 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))))) * x4 + (( (2/1) * x7 + (x10 + (x21 + (x22 + (x23 + (x24 + (x25 + (x26 + (x27 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))))) * x5 + (( (-2/1) * x7 + ( - x11 + ( - x43 + ( - x44 + ( - x45 + ( - x48)))))) * x6))) * x2 + (((x7 + (x8 + (x9 + (x13 + (x16 + (x17 + (x18 + (x19 + (x20 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))))) * x4 + (( - x7 + ( - x8 + ( - x10 + ( - x40 + ( - x41 + ( - x42 + ( - x45 + ( - x48)))))))) * x5 + ((x7 + (x8 + (x11 + (x29 + (x30 + (x31 + (x32 + (x33 + (x34 + (x35 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))))) * x6))) * x3 + ((( - x7 + ( (-2/1) * x9 + ( - x13 + ( - x17 + ( - x18 + ( - x19 + ( - x20 + ( - x28 + ( - x36 + ( - x39 + ( - x42 + ( - x45 + ( - x48))))))))))))) * x4 + ((x7 + (x9 + (x10 + (x13 + (x18 + (x19 + (x20 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))) * x5 + ((x7 + (x9 + (x11 + (x13 + (x19 + (x20 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48)))))))))))) * x6))) * x4)))) * x1 + (((( (-2/1) * x7 + ( - x10 + ( - x21 + ( - x22 + ( - x23 + ( - x24 + ( - x25 + ( - x26 + ( - x27 + ( - x28 + ( - x36 + ( - x39 + ( - x42 + ( - x45 + ( - x48))))))))))))))) * x5) * x2 + ((( - x7 + ( - x8 + ( - x9 + ( - x37 + ( - x38 + ( - x39 + ( - x42 + ( - x45 + ( - x48))))))))) * x4 + ((x7 + (x8 + (x10 + (x21 + (x23 + (x24 + (x25 + (x26 + (x27 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))))) * x5 + ((x7 + (x8 + (x11 + (x29 + (x30 + (x31 + (x32 + (x33 + (x34 + (x35 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))))) * x6))) * x3 + (((x7 + (x9 + (x10 + (x21 + (x24 + (x25 + (x26 + (x27 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48)))))))))))))) * x5) * x4 + ((( - x7 + ( (-2/1) * x10 + ( - x21 + ( - x25 + ( - x26 + ( - x27 + ( - x28 + ( - x36 + ( - x39 + ( - x42 + ( - x45 + ( - x48)))))))))))) * x5 + ((x7 + (x10 + (x11 + (x21 + (x26 + (x27 + (x28 + (x36 + (x39 + (x42 + (x45 + (x48)))))))))))) * x6)) * x5)))) * x2 + (((( (-2/1) * x8 + ( - x11 + ( - x29 + ( - x31 + ( - x32 + ( - x33 + ( - x34 + ( - x35 + ( - x36 + ( - x39 + ( - x42 + ( - x45 + ( - x48))))))))))))) * x6) * x3 + (((x8 + (x9 + (x11 + (x29 + (x32 + (x33 + (x34 + (x35 + (x36 + (x39 + (x42 + (x45 + (x48))))))))))))) * x6) * x4 + (((x8 + (x10 + (x11 + (x29 + (x33 + (x34 + (x35 + (x36 + (x39 + (x42 + (x45 + (x48)))))))))))) * x6) * x5 + (( - x8 + ( (-2/1) * x11 + ( - x29 + ( - x34 + ( - x35 + ( - x36 + ( - x39 + ( - x42 + ( - x45 + ( - x48)))))))))) * pow(x6,2))))) * x3 + (((( - x9 + ( - x10 + ( - x11 + ( - x46 + ( - x47 + ( - x48)))))) * x6) * x5) * x4)));
//	eq = eq.expand();
////	cout << "equation :: " << eq << endl;
//	lst params;
//	lst vars;
//	vars = x1,x2,x3,x4,x5,x6;
//	params = 	x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,
//				x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,
//				x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,
//				x40,x41,x42,x43,x44,x45,x46,x47,x48;
//	// Espace des paramètres
//
//
//	vector<double> pAi (params.nops(),0);
//	vector< vector<double> > pA (2*params.nops(),pAi);
//	vector<double> pb (2*params.nops(),0);
//
//	// left constraint .... <= g(x) --> -g(x) <= -....
//
//	for(int i=0;i<2*params.nops();i++)
//	{
//		if(i<(params.nops()))
//			pA[i][i] = -1;
//		else
//			pA[i][i-(params.nops())] = 1;
//
//		pb[i] = 1;
//	}
//
////	cout << "parameters lin sys\n";
//	LinearSystem *parameters = new LinearSystem(pA,pb);
//	LinearSystemSet *paramSet = new LinearSystemSet(parameters);
//
//	int dim_global = 6; // Dimension of the complete system (not just this equation)
//	// Constructeur Equation
//
//	// Calcul with parallelotope description (here a [0;1] box)
//	// DEFINITION Of the parallelotopes
//	possymbol q1("q1"), q2("q2"), q3("q3"), q4("q4"), q5("q5"), q6("q6"), q7("q7");
//	possymbol a1("a1"), a2("a2"), a3("a3"), a4("a4"), a5("q5"), a6("q6"), a7("q7");
//	possymbol b1("b1"), b2("b2"), b3("b3"), b4("b4"), b5("q5"), b6("q6"), b7("q7");
//	lst qs,as,bs; // Symbolic variables for linear transformation
//	qs = q1,q2,q3,q4,q5,q6;
//	as = a1,a2,a3,a4,a5,a6;
//	bs = b1,b2,b3,b4,b5,b6;
//	vector<lst> set_vars;
//	set_vars.push_back(qs);
//	set_vars.push_back(as);
//	set_vars.push_back(bs);
//	doubleVector row = doubleVector(6,0);
//	vector<doubleVector> u= vector<doubleVector>(6,row);
//	u[0][0] = 1;
//	u[1][1] = 1;
//	u[2][2] = 1;
//	u[3][3] = 1;
//	u[4][4] = 1;
//	u[5][5] = 1;
//	Parallelotope* para = new Parallelotope(set_vars,u);
//
//	doubleVector bV = {4,4,4,4,4,4};
//	doubleVector ampl = {2.36,2.36,2.36,2.36,2.36,2.36};
//
//	poly_values p;
//	p.base_vertex = bV;
//	p.versors = u;
//	p.lenghts = ampl;
//	para->add_polyValues(p);
//
//	// Computing f°g
//	lst generator_function = para->getGeneratorFunction();
//
//	lst sub;
//	ex fog;
//	for(int i=0; i<(signed)generator_function.nops(); i++){
//		sub.append(vars[i] == generator_function[i]);
//	}
//	fog = eq.subs(sub,subs_options::no_pattern);
//
//	Equation* eq_p = new Equation(	dim_global, as, params, fog, EXPLICIT_SYM, qs, bs);
//
//	tStart = clock();
//	pair<double,double>  res_explicit = eq_p->optimize(para,paramSet);
//	cout << "Time method explicit Parallelotopes: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
//	cout << "Result from optimization of: \n" << eq << endl;
//	cout << "min = " << -1*res_explicit.first*err << endl;
//	cout << "max = " << res_explicit.second*err << endl;
//

}





