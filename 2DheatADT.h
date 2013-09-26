/*******************************************************************
 * file:    2DheadADT.h
 * version: 1.0
 * project: Comp7751 Assignment Project
 * author:  Song Gao
 * date:    April 25, 2003
 * description: Program to solve the 2D Diffusion (heat) equation by 
				using the Crank-Nicolson-ADT scheme. 
				Under Neumann boundary conditions:
				u_x(0, y, t)=0,		u_x(L, y, t)=0,
				u_x(x, 0, t)=0,		u_y(x, L, t)=0.
	
				Initial condition is:
				u(x, y, t=0) = f(x,y)
 * 	           
 ******************************************************************/ 

#ifndef HEATADT_CLASS
#define HEATADT_CLASS

#include "DMatrix.h"
#include "Image.h"


class heatADT 
{
  private:
	int nx;
	int ny;
	double dx;
	double dy;
//	Matrix u;
  public:
	heatADT();
	~heatADT(){};
	void init(Matrix&, Image&); 
	Matrix getrhsB(int, double, Matrix&);
	Matrix StencilB(double&);
	Matrix getrhsQ(int, double, Matrix&);
	Matrix StencilQ(double);
	
	void tridge( Matrix, Matrix, Matrix&);
	Matrix solver(double kapp, double dt, Matrix &u);
	Matrix solveDiffu(double kapp, double dt, Matrix& u);
	Matrix FTCSsolver(double mu, double h, double dt, Matrix& u,
					  Image&, double&, double&);
};

#endif  
