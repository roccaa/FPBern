/*
 * Parallelotope.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: dreossi
 */

#include "Parallelotope.h"

Parallelotope::Parallelotope(vector<lst> vars, vector< vector<double> > u) {

	if(vars.size() != 3){
		cout<<"Parallelotope::Parallelotope : vars must contain 3 collections of variable names (q,alpha,beta)";
		exit (EXIT_FAILURE);
	}
	this->re_used = false;
	this->vars.push_back(vars[0]);
	this->vars.push_back(vars[1]);
	this->vars.push_back(vars[2]);

	// get the dimension of the parallelotope
	this->dim = vars[0].nops();

	// and store its variable names
	for(int i=0; i<3; i++){
		if((signed)vars[i].nops() != this->dim){
			cout<<"Parallelotope::Parallelotope : vars["<<i<<"] must have "<<this->dim<<" variables";
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

	// create the generation function accumulating the versor values
	for(int i=0; i<this->dim; i++){

		// check dimension of versors
		if((signed)u[i].size() != this->dim){
			cout<<"Parallelotope::Parallelotope : dim and ui dimensions must agree";
			exit (EXIT_FAILURE);
		}

		// Normalize the versors and generatre the generator function
		double norm = euclidNorm(u[i]);
		vector< double > norm_versor;

		for(int j=0; j<this->dim; j++){
			norm_versor.push_back(u[i][j]/norm);
			double normVersor_c = u[i][j]/norm;
			this->generator_function[j] = this->generator_function[j] + alpha[i]*beta[i]*normVersor_c;
		}
		// store the i-th versor
		this->u.push_back(norm_versor);
	}

	// Initialize the template matrix
	vector< double > base_vertex (this->dim,0);
	vector<double> lenghts (this->dim,1);

	//TODO TEMPORAIRE
	//LinearSystem *LS = this->gen2const(base_vertex,lenghts);
	//this->template_matrix = LS->getA();
	this->linSystem = NULL;

	this->type = PARALLELOTOPE;
	this->overBox = NULL;
}

Parallelotope::Parallelotope(vector<lst> vars, LinearSystem *constr) {

	//cout << "Parallelotope constructor from linear system\n";
	this->overBox = NULL;
	if(vars.size() != 3){
		cout << "vars.size() = " << vars.size() << endl;
		cout<<"Parallelotope::Parallelotope : vars must contain 3 collections of variable names (q,alpha,beta)";
		exit (EXIT_FAILURE);
	}

	this->vars.push_back(vars[0]);
	this->vars.push_back(vars[1]);
	this->vars.push_back(vars[2]);

	// get the dimension of the parallelotope
	this->dim = vars[0].nops();
	// and store its variable names
	for(int i=0; i<3; i++){
		if((signed)vars[i].nops() != this->dim){
			cout<<"Parallelotope::Parallelotope : vars["<<i<<"] must have "<<this->dim<<" variables";
			exit (EXIT_FAILURE);
		}
	}

	// initialize generator function
	for(int i=0; i<this->dim; i++){
		this->generator_function.append(this->vars[0][i]);
	}

	// convert the linear system to vectors
	vector< vector<double> > Lambda = constr->getA();
	vector<double> d = constr->getb();
	vector< vector<double> > vertices;

	// find base vertex
	//build the linear system
	ex q = this->vars[0];
	lst LS;
//	cout << "Lambda = \n";
//	for(int i=0; i<Lambda.size(); i++){
//		for(int j=0; j<Lambda[i].size(); j++){
//			cout<<Lambda[i][j]<<" ";
//		}
//		cout<<"\n";
//	}
//	cout << "d = " << d;

	for(int i=0; i<this->dim; i++){
		ex eq = 0;
		for(int j=0; j<(signed)Lambda[i].size(); j++){
			eq = eq + Lambda[i][j]*q[j];
		}
		eq = eq == d[i];
		LS.append(eq);
	}
//	cout << "Linear system to solve to find base vertex :: \n";
//	cout<<LS;
//	cout << "\n";

	ex solLS = lsolve(LS,q);
//	cout << "solution is:\n";
//	cout << solLS << endl;
	if( solLS.nops() == 0 ){	// the template is singular
		cout<<"singluar parallelotope\n";
		constr->print();
		cout<<LS;
		assert(0);
		return;
	}

	vertices.push_back( lst2vec((ex)q.subs(solLS)) ); // store the base_vertex
//	cout << "base_vertex is :: " << vertices[0];
	doubleVector bl = vertices[0];
	doubleVector bu = vertices[0];
	int n = 0;

	// Compute the vertices v
	for(int k=0; k<this->dim; k++){
		ex a = this->vars[1];
		lst LS;

		for(int i=0; i<this->dim; i++){
			ex eq = 0;
			for(int j=0; j<(signed)Lambda[i].size(); j++){
				eq = eq + Lambda[i][j]*a[j];
			}
			if(i != k){
				eq = eq == d[i];
			}else{
				eq = eq == -d[i+this->dim];
			}
			LS.append(eq);
		}
//		cout << "##### Next vertex computation \n";
//		cout << "solLS = \n";
//		cout << solLS << endl;
		ex solLS = lsolve(LS,a);

		vertices.push_back( lst2vec(a.subs(solLS)) ); // store the i-th vertex
//		cout << "new vertex is  = " << vertices[k+1] << endl;
		n++;
		for(int j=0;j<this->dim;j++)
		{
			if(vertices[n][j]<bl[j])
				bl[j] = vertices[n][j];
			if(vertices[n][j]>bu[j])
				bu[j] = vertices[n][j];
		}
	}

	this->actual_base_vertex = vertices[0];

	// Compute the generators
//	cout << "Compute the generators\n";
	vector< vector<double> > g;
	for(int i=0; i<this->dim; i++){
		vector<double> gi;
		for(int j=0; j<this->dim; j++){
			gi.push_back( vertices[i+1][j] - vertices[0][j] );
		}
//		cout << "g" << i <<" = " << gi;
		g.push_back(gi);
	}

	// Compute the generators lengths
//	cout << "length = ";
	vector< double > lengths;
	for(int i=0; i<this->dim; i++){
		lengths.push_back(euclidNorm(g[i]));
//		cout<<lengths[i]<<" ";
	}
//	cout<<"\n";

	this->actual_lenghts = lengths;

//	cout<<"g\n";
	// Find the versors
	vector< vector< double > > versors;
	for(int i=0; i<this->dim; i++){
		vector<double> versori;
		for(int j=0; j<this->dim; j++){
//			cout<<g[i][j]<<" ";
			if(lengths[i] != 0){
				versori.push_back( floor((g[i][j]/lengths[i]) * 100000000000.0f) / 100000000000.0f );
			}else{
				versori.push_back( floor(g[i][j] * 100000000000.0f) / 100000000000.0f );
			}
		}
//		cout<<"\n";
		versors.push_back(versori);
	}

	this->u = versors;
//	cout << "u = \n";
//	for(int i=0; i<this->u.size(); i++){
//		for(int j=0; j<this->u[i].size(); j++){
//			cout<<this->u[i][j]<<" ";
//		}
//		cout<<"\n";
//	}
//	cout << "lengths = " << lengths;
//	cout << "base vertex = " << vertices[0];

	// create the generation function accumulating the versor values
	lst alpha = this->vars[1];
	lst beta = this->vars[2];
//	cout << "alpha = " << alpha << endl;
//	cout << "beta = " << beta << endl;
	for(int i=0; i<this->dim; i++){

		// check dimension of versors
		if((signed)u[i].size() != this->dim){
			cout<<"Parallelotope::Parallelotope : dim and ui dimensions must agree";
			exit (EXIT_FAILURE);
		}

		// Generatre the generator function
		//double norm = euclidNorm(u[i]);
		//vector< double > norm_versor;
//		cout << "u[i] = " << u[i];
		for(int j=0; j<this->dim; j++){
			//norm_versor.push_back(u[i][j]/norm);
			this->generator_function[j] = this->generator_function[j] + alpha[i]*beta[i]*this->u[i][j];
//			cout << "this->generator_function[j] = " << this->generator_function[j] << endl;
		}
	}
//	cout << "generator_function = " << this->generator_function << endl;
//	assert(0);
	this->template_matrix = constr->getA();
	this->type = PARALLELOTOPE;
	this->linSystem = constr;

	doubleVector pointe = vertices[0];
	for(int i=0;i<versors.size();i++)
	{
		pointe = pointe + lengths[i]*versors[i];
	}
	vertices.push_back(pointe);

//	cout << "Vertices:\n";
//	for(int i=0;i<vertices.size();i++)
//	{
//		cout << "v" <<i<<" = " << vertices[i];
//	}

	this->vertices_list = vertices;

	for(int j=0;j<this->dim;j++)
	{
		if(pointe[j]<bl[j])
			bl[j] = pointe[j];
		if(pointe[j]>bu[j])
			bu[j] = pointe[j];
	}

	this->num_genDescription.base_vertex = vertices[0];
	this->num_genDescription.lenghts = lengths;
	this->num_genDescription.versors = this->u;
	if(this->overBox != NULL)
			delete(this->overBox);
	this->overBox = new Box(bl,bu);

}



vector_matrix Parallelotope::generate_template(vector_matrix normals)
{
	vector_matrix res;
	for(int i=0;i<normals.size();i++)
	{
		res.push_back(normals[i]);
	}
	for(int i=0;i<normals.size();i++)
	{
		res.push_back(-1*normals[i]);
	}
	return res;

}


void Parallelotope::recompute_gen_func(vector_matrix u)
{
	this->u.clear();
	lst q = this->vars[0];
	lst alpha = this->vars[1];
	lst beta = this->vars[2];
	// initialize generator function
	lst empty;
	this->generator_function = empty;
	for(int i=0; i<this->dim; i++){
		this->generator_function.append(q[i]);
	}

	// create the generation function accumulating the versor values
	for(int i=0; i<this->dim; i++){

		// check dimension of versors
		if((signed)u[i].size() != this->dim){
			cout<<"Parallelotope::Parallelotope : dim and ui dimensions must agree";
			exit (EXIT_FAILURE);
		}

		// Normalize the versors and generatre the generator function
		double norm = euclidNorm(u[i]);
		vector< double > norm_versor;

		for(int j=0; j<this->dim; j++){
			norm_versor.push_back(u[i][j]/norm);
			double normVersor_c = u[i][j]/norm;
			//this->generator_function[j] = this->generator_function[j] + alpha[i]*beta[i]*u[i][j];
			this->generator_function[j] = this->generator_function[j] + alpha[i]*beta[i]*normVersor_c;
		}
		// store the i-th versor
		//this->u.push_back(u[i]);
		this->u.push_back(norm_versor);
	}
}


void  Parallelotope::add_polyValues(poly_values p)
{
	this->num_genDescription = p;
	//cout << p.base_vertex;
	if(this->overBox==NULL)
		this->compute_overBox();
//	clear_BCinfo();

	//	cout << "##new overBox:\n" << this->overBox->getBl() << this->overBox->getBu();

	//	cout << "there is " << this->vertices_list.size() << " vertices\n";
//	cout << "the vertices are: \n";
//	for(int i =0;i<this->vertices_list.size();i++)
//	{
//		cout << vertices_list[i] << ";";
//	}

}

void  Parallelotope::add_polyValues(poly_values p, bool do_overBox)
{
	cout << "add_poly_value V2\n";
	this->num_genDescription = p;
	if((this->overBox==NULL) || (do_overBox==true))
		this->compute_overBox();
//	cout << "BC_info size = " << this->info_BernsteinCoeff.size();
	clear_BCinfo();
	cout << "##new overBox:\n" << this->overBox->getBl() << this->overBox->getBu();
	cout << "ther vertices are: \n";
	for(int i =0;i<this->vertices_list.size();i++)
	{
		cout << vertices_list[i] << ";";
	}
}

Box * Parallelotope::compute_overBox()
{
//	cout << "compute over box\n";
//	cout << "base vertex = " << this->num_genDescription.base_vertex;
//	cout << "length = " << this->num_genDescription.lenghts;
	LinearSystem *linsyst = gen2const(this->num_genDescription.base_vertex,this->num_genDescription.lenghts);
//	cout << "system is \n";
//	linsyst->print();
	vector< vector<double> > Lambda = linsyst->getA();
	vector<double> d = linsyst->getb();
	vector< vector<double> > vertices;

	ex q = this->vars[0];
	lst LS;

	for(int i=0; i<this->dim; i++){
		ex eq = 0;
		for(int j=0; j<(signed)Lambda[i].size(); j++){
			eq = eq + Lambda[i][j]*q[j];
		}
		eq = eq == d[i];
//		cout << "eq = " << eq << endl;
		LS.append(eq);
	}

	ex solLS = lsolve(LS,q);
	//cout << "solution symbolic = " << solLS << endl;
	vertices.push_back( lst2vec(q.subs(solLS,subs_options::no_pattern)) ); // store the base_vertex
	doubleVector bl = vertices[0];
	doubleVector bu = vertices[0];
	int n=0;

	// Compute the vertices v
	for(int k=0; k<this->dim; k++){
		ex a = this->vars[1];
		lst LS;

		for(int i=0; i<this->dim; i++){
			ex eq = 0;
			for(int j=0; j<(signed)Lambda[i].size(); j++){
				eq = eq + Lambda[i][j]*a[j];
			}
			if(i != k){
				eq = eq == d[i];
			}else{
				eq = eq == -d[i+this->dim];
			}
			LS.append(eq);
		}

		ex solLS = lsolve(LS,a);
		vertices.push_back( lst2vec(a.subs(solLS,subs_options::no_pattern)) ); // store the i-th vertex
		n++;
		for(int j=0;j<this->dim;j++)
		{
			if(vertices[n][j]<bl[j])
				bl[j] = vertices[n][j];
			if(vertices[n][j]>bu[j])
				bu[j] = vertices[n][j];
		}
//		cout << "vertices["<<n<<"] = " << vertices[n];
	}
	doubleVector pointe = vertices[0];
	for(int i=0;i<this->num_genDescription.versors.size();i++)
	{
		pointe = pointe + this->num_genDescription.lenghts[i]*this->num_genDescription.versors[i];
	}
//	cout << "pointe vertex = " << pointe ;
	vertices.push_back(pointe);
	this->vertices_list = vertices;
	for(int j=0;j<this->dim;j++)
	{
		if(pointe[j]<bl[j])
			bl[j] = pointe[j];
		if(pointe[j]>bu[j])
			bu[j] = pointe[j];
	}
	//assert(this->overBox==NULL);
	if(this->overBox != NULL)
		delete(this->overBox);
	this->overBox = new Box(bl,bu);
	// TODO TEMPORAIRE
	//linsyst->print();
	if(this->linSystem == NULL)
	{
		this->linSystem = linsyst;
		this->template_matrix = linsyst->getA();
	}
	else
		delete(linsyst);
	return this->overBox;

}

// Convert from generator to constraints representation
// q : numeric base vertex
// beta : numeric lenghts
LinearSystem* Parallelotope::gen2const(vector<double> q, vector<double> beta){

	if((signed)q.size() != this->dim){
		cout<<"Parallelotope::gen2const : q must have dimension "<<this->dim;
		exit (EXIT_FAILURE);
	}


	vector< vector<double> > hps;	// hyperplane equations

//	cout << "First set of hyperplanes\n";
	// first set of hyperplanes
	for(int i=0; i<this->dim; i++){
//		cout << "i = " << i << endl;
		vector< vector<double> > pts;
		pts.push_back(q);							// add base vertex
//		cout << "added pt = " << q ;
		for(int j=0; j<this->dim; j++){ 			// for all the generators u
			if( i != j ){
				vector<double> p;
				for(int k =0; k<this->dim; k++){	// coordinate of the point
					p.push_back(q[k] + this->u[j][k]*beta[j]);
				}
				pts.push_back(p);
//				cout << "added pt = " << p ;
			}
		}
		hps.push_back(hyperplaneThroughPts(pts));
	}

//	cout << "\nSecond set of hyperplanes\n";
	//second set of hyperplanes
	for(int i=0; i<this->dim; i++){

			vector< vector<double> > pts;

			vector<double> qt; // traslate q
			for(int j=0; j<this->dim; j++){
				qt.push_back(q[j] + beta[i]*this->u[i][j]);
			}
			pts.push_back(qt);							// add base vertex
//			cout << "added pt = " << qt ;
			for(int j=0; j<this->dim; j++){ 			// for all the generators u
				if( i != j ){
					vector<double> p;
					for(int k =0; k<this->dim; k++){	// coordinate of the point
						p.push_back( (q[k] + this->u[j][k]*beta[j]) + beta[i]*this->u[i][k]);
					}
					pts.push_back(p);
//					cout << "added pt = " << p ;
				}
			}

			hps.push_back(hyperplaneThroughPts(pts));

	}

	vector< vector<double> > Lambda;
	vector< double > d (this->dim*2,0);

	// initialize template Lambda
	vector<double> Lambda_i (this->dim,0);
	for(int i=0; i<(this->dim)*2; i++){
		Lambda.push_back(Lambda_i);
	}

	for(int i=0; i<this->dim; i++){

		// find hyperplane with smaller direction
		double b1 = -hps[i][hps[i].size()-1];
		double b2 = -hps[i+this->dim][hps[i].size()-1];

		if(b1 > b2){
			d[i] = b1;
			d[i+this->dim] = -b2;
			for(int j=0; j<this->dim; j++){
				Lambda[i][j] = hps[i][j];
				Lambda[i+this->dim][j] = -hps[i+this->dim][j];
			}
		}else{
			d[i] = -b1;
			d[i+this->dim] = b2;
			for(int j=0; j<this->dim; j++){
				Lambda[i][j] = -hps[i][j];
				Lambda[i+this->dim][j] = hps[i+this->dim][j];
			}
		}
	}

	LinearSystem *LS = new LinearSystem(Lambda,d);

	return LS;
}


// Compute the equation of the hyperplane passing through the points in pts
// the equation is res[0]*x_0 + .. + res[n]*x_n + res[n+1] = 0
// pts : matrix containing the points
vector<double> Parallelotope::hyperplaneThroughPts(vector< vector<double> > pts){

	if(pts.size() != pts[0].size()){
		cout<<"Parallelotope::hyperplaneThroughPts : pts must contain n n-dimensional points";
		exit (EXIT_FAILURE);
	}

	vector< vector<double> > A;
	// build the linear system Ax = 0
	for(int i=1; i<(signed)pts.size(); i++){
		vector<double> Ai;
		for(int j=0; j < this->dim; j++){
			Ai.push_back(pts[0][j] - pts[i][j]);
		}
//		cout << "base = " << pts[0];
//		cout << "pts[" << i << "] = " << pts[i];
//		cout << "base - pts[" << i << "] = " << Ai;
		A.push_back(Ai);
	}

	// Build the linear system to find the normal vector
	lst LS;
	lst a = this->vars[1];
//	cout << "(Ax = 0) LS = \n";
	for(int i=0; i<this->dim-1; i++){

		ex eq = 0;
		lst sub;

		for(int j=0; j<this->dim; j++){
			eq = eq + A[i][j]*a[j];
		}

		eq = eq == 0;
		LS.append(eq);
//		cout << eq << endl;
	}

	// Solve the linear system
//	cout << "Solve the linear system\n";
	ex solLS = lsolve(LS,a);
//	cout << solLS << endl;
//	cout << "a = " << a << endl;
//	cout << "done !\n";

	// Build the linear inequality
	lst x = this->vars[0];
	ex eqr = 0;
	vector<double> lambda (this->dim+1,0);

	ex sub;
	int sub_idx;
	// search for the tautology
	for(int i=0; i<(signed)solLS.nops(); i++){
		if(solLS[i].is_equal(a[i] == a[i])){
//			cout << "Tautology\n";
//			cout << a[i] << " = " << 1 << endl;
			sub = a[i] == 1;
			sub_idx = i;
		}
	}

	lambda[sub_idx] = 1;
	for(int i=0; i<this->dim; i++){
		if(i != sub_idx){
//			cout << "subs result = " << solLS[i].subs(sub,subs_options::no_pattern) << endl;
//			cout << "a[i] = " << a[i] << endl;
//			cout << "a[i].subs() = " << a[i].subs(solLS[i].subs(sub,subs_options::no_pattern)) << endl;
			lambda[i] = ex_to<numeric>(evalf(a[i].subs(solLS[i].subs(sub,subs_options::no_pattern)))).to_double();
		}
		eqr = eqr + a[i].subs(a[i] == lambda[i],subs_options::no_pattern)*pts[0][i];
	}

	lambda[this->dim] = -ex_to<numeric>(evalf(eqr)).to_double();

//	cout << "Lambda = " << lambda << endl;
	return lambda;
}


// Convert the constrain representation to the generator one
// constr : LinearSystem
poly_values Parallelotope::const2gen(LinearSystem *constr){

	//cout << "Parallelotope::const2gen\n";
	vector< vector<double> > Lambda = constr->getA();
	vector<double> d = constr->getb();
	vector< vector<double> > vertices;

	//cout << "Linear system:\n";
	//constr->print();

	// find base vertex
	//build the linear system
	ex q = this->vars[0];
	//cout << "q = " << q << endl;
	lst LS;

	for(int i=0; i<this->dim; i++){
		ex eq = 0;
		for(int j=0; j<(signed)Lambda[i].size(); j++){
			eq = eq + Lambda[i][j]*q[j];
		}
		eq = eq == d[i];
		LS.append(eq);
	}

	ex solLS = lsolve(LS,q);
	if( solLS.nops() == 0 ){	// the template is singular
		cout<<"singluar parallelotope\n";
		constr->print();
		cout<<LS;
		assert(0);
	}

	vertices.push_back( lst2vec(q.subs(solLS,subs_options::no_pattern)) ); // store the base_vertex

	doubleVector bl = vertices[0];
	doubleVector bu = vertices[0];
	int n=0;

	// Compute the vertices v
//	cout << "Compute the vertices\n";
	for(int k=0; k<this->dim; k++){
		ex a = this->vars[1];
		lst LS;

		for(int i=0; i<this->dim; i++){
			ex eq = 0;
			for(int j=0; j<(signed)Lambda[i].size(); j++){
				eq = eq + Lambda[i][j]*a[j];
			}
			if(i != k){
				eq = eq == d[i];
			}else{
				eq = eq == -d[i+this->dim];
			}
			LS.append(eq);
		}

		ex solLS = lsolve(LS,a);
		vertices.push_back( lst2vec(a.subs(solLS,subs_options::no_pattern)) ); // store the i-th vertex
		n++;
		for(int j=0;j<this->dim;j++)
		{
			if(vertices[n][j]<bl[j])
				bl[j] = vertices[n][j];
			if(vertices[n][j]>bu[j])
				bu[j] = vertices[n][j];
		}
	}


	// Compute the generators
//	cout << "Compute the generators \n";
	vector< vector<double> > g;
	for(int i=0; i<this->dim; i++){
		vector<double> gi;
		for(int j=0; j<this->dim; j++){
			gi.push_back( vertices[i+1][j] - vertices[0][j] );
		}
		g.push_back(gi);
	}

	// Compute the generators lengths
//	cout << "Compute the generators lengths \n";
	vector< double > lengths;
	for(int i=0; i<this->dim; i++){
		//cout << "generator = " << g[i];
		double lgt = euclidNorm(g[i]);
		lengths.push_back(lgt);
	}

	// Find the versors
//	cout << "Find the versors \n";
	vector< vector< double > > versors;
	for(int i=0; i<this->dim; i++){
		vector<double> versori;
		// TODO here its different from the LS constructor ...
		// Do not take into account the case where lengthi = 0
		for(int j=0; j<this->dim; j++){
			versori.push_back( g[i][j]/lengths[i] );
		}
		versors.push_back(versori);
	}

	// Return the conversion
	poly_values result;
	result.base_vertex = vertices[0];
	result.lenghts = lengths;
	result.versors = versors;
	doubleVector pointe = vertices[0];

	for(int i=0;i<result.versors.size();i++)
	{
		pointe = pointe + lengths[i]*versors[i];
	}
	vertices.push_back(pointe);

	this->vertices_list = vertices;
	for(int j=0;j<this->dim;j++)
	{
		if(pointe[j]<bl[j])
			bl[j] = pointe[j];
		if(pointe[j]>bu[j])
			bu[j] = pointe[j];
	}
	if(this->overBox != NULL)
			delete(this->overBox);
	this->overBox = new Box(bl,bu);

	return result;
}

// Convert a symbolic list to a vector
vector<double> Parallelotope::lst2vec(ex list){

	vector< double > res;

	for(int i=0; i<(signed)list.nops(); i++)
		res.push_back( ex_to<numeric>(evalf( list[i] )).to_double() );

	return res;
}

// Compute the Euclidean norm of the vector v
double Parallelotope::euclidNorm(vector<double> v){

	double norm = 0;

	for(int i=0; i<(signed)v.size(); i++){
		norm = norm + (v[i]*v[i]);
	}

	return sqrt(norm);
}

// TODO DELETE
Box* Parallelotope::compute_reachedBox(void){
//	assert(this->info_BernsteinCoeff.empty()==false);
//	doubleVector bl = doubleVector(this->info_BernsteinCoeff.size());
//	doubleVector bu = doubleVector(this->info_BernsteinCoeff.size());
////	assert(this->info_BernsteinCoeff.size() == dim); // Meme nombre dynamiques et de variables ...
//	for(int i=0; i<this->info_BernsteinCoeff.size(); i++){
//		bl[i] = this->info_BernsteinCoeff[i]->min;
//		bu[i] = this->info_BernsteinCoeff[i]->max;
//	}
//	return new Box(bl,bu);

}
/*
void Parallelotope::clear_BCinfo(void){
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

void Parallelotope::clear_BCinfo(void)
{

}

Parallelotope::~Parallelotope() {
	// TODO Auto-generated destructor stub
	delete(this->overBox);

	clear_BCinfo(); // TODO  remove this struct as soon as possible

	if(this->linSystem!=NULL)
		delete(this->linSystem);
}

