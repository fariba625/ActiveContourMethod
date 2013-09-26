

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
