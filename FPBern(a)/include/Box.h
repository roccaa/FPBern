/*
 * Box.h
 *
 *  Created on: Nov 3, 2014
 *      Author: dreossi
 */

#ifndef BOX_H_
#define BOX_H_

#include "Polyhedron.h"

class Box : public Polyhedron{
private:
	vector<double> l;
	vector<double> u;
	int LocId;
	int motherLocId;
	// For each equation there is one Berstein Coeff structure



public:
//	vector<BernCoeff_for_Box *> info_BernsteinCoeff;
	int dim;
	Box(void){};
	Box(vector<double> l, vector<double> u);
	Box(LinearSystem* linSystem);
	Box(vector<lst> vars);
	LinearSystem* bounds2const(vector<double> l, vector<double>u);
	vector< vector<double> > const2bounds(LinearSystem* linSystem);
	LinearSystem* gen2const(vector<double> q, vector<double> beta);
	poly_values const2gen(LinearSystem *constr);
	vector<double> getBu(void){return this->u;};
	vector<double> getBl(void){return this->l;};
	double compute_directionSize(int direction){return std::abs(u[direction]-l[direction]);};
	void setLocationId(int id){this->LocId = id;};
	int getLocationId(void){return this->LocId;};
	void setMotherLocationId(int id){this->motherLocId = id;};
	int getMotherLocationId(void){return this->motherLocId;};
	double get_fullBoxSize();
	int compute_MaxSizeDirection();
/*	void add_BernsteinInfo(BernCoeff_for_Box * info)
	{
		info_BernsteinCoeff.push_back(info);
	};*/
//	vector<BernCoeff_for_Box *> *  get_bernsteinInfo(){return &info_BernsteinCoeff;}
	Box *compute_reachedBox(void);
	void reset_lin_syst(void){this->linSystem = NULL;};
//	void clear_BCinfo(void);
	virtual ~Box();

};

#endif /* BOX_H_ */
