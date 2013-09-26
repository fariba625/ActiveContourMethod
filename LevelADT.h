/*******************************************************************
 * file:    LevelADT.h
 * version: 1.0
 * project:  Project
 * author:  Song Gao
 * date:    May 01, 2003
 * description: Program to solve the Level set equation by 
				using the Crank-Nicolson-ADT scheme. 
				Under Neumann boundary conditions:
				u_x(0, y, t)=0,		u_x(L, y, t)=0,
				u_x(x, 0, t)=0,		u_y(x, L, t)=0.
	
				Initial condition is:
				u(x, y, t=0) = f(x,y)
 * 	           
 ******************************************************************/ 

#ifndef LEVELADT_CLASS
#define LEVELADT_CLASS

#include "DMatrix.h"
#include "Image.h"


class LevelADT 
{
  private:
	int nx;
	int ny;
	double dx;
	double dy;
//	Matrix u;
  public:
	LevelADT();
	LevelADT(int,int);
	~LevelADT(){};
//	void init(Matrix&); 
	Matrix getrhsB(int k, Matrix &ry1, Matrix &ry2, Matrix &DT, Matrix &u);
	Matrix StencilB(int j, Matrix &rx1, Matrix &rx2);
	Matrix getrhsQ(int j, Matrix &rx1, Matrix &rx2, Matrix &DT, Matrix &u);
	Matrix StencilQ(int j, Matrix &ry1, Matrix &ry2);
	
	void tridge( Matrix, Matrix, Matrix&);
	Matrix solver(double mu, double dt, Matrix& u,Image &U0,double &U1, double &U2);

    void FTCSsolver(double mu,double dt, Matrix& u,Image&, double&, double&);
	void phi_1solver(double, double, Matrix&, Image&, double, double, 
					   double, double, Matrix&);
	void phi_2solver(double, double, Matrix&, Image&, double, double, 
					   double, double, Matrix&);
	Matrix ADTsolver1(double, double, Matrix&, Image&, double, double, 
					   double, double, Matrix&);
	Matrix ADTsolver2(double, double, Matrix&, Image&, double, double, 
					   double, double, Matrix&);

};

#endif  
