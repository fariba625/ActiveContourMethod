/*******************************************************************
 * file:    2DheadADT.cpp 
 * version: 1.0
 * project: Comp7751 Assignment Project
 * author:  Song Gao
 * date:    April 25, 2003
 * description: Function implementations in class heatADT
 * 	           
 ******************************************************************/ 

#include <iostream>
#include <fstream>
#include <assert.h>  
#include <math.h>
#include "DMatrix.h"

#include "2DheatADT.h"

#define PI 3.1415926

         
// default constructor
heatADT::heatADT()
{    
	Image temp;
	nx = temp.getRow(); //256; 
	ny = temp.getCol(); //256; 
	dx = 1.;
	dy =1.;
}

// Initial conditions 
void heatADT::init(Matrix &u, Image &u0)
{
	int j,k;

	for(j=1; j<=nx; j++)
	for(k=1; k<=ny; k++)
	{ // a circle with Radius = 80
	//	u(j, k) = sqrt(pow((j - (double)nx/2.0),2)+ 
	//			  pow((k - (double)ny/2.0),2))-100.; //80
		u(j,k) = u0(j,k);
	}
}


// this function return the right hand side of the first equation
// here u is u^{n+1/2}
Matrix heatADT::getrhsB(int k, double ry, Matrix &u)
{
	int j,  kk=1;
	Matrix rhSide(nx,1);
	
		for(j=1; j<nx-1; j++) 
		{
			if (k == 1)
			{
				rhSide(kk,1) = ry*u(j, k+1) + (1.- ry)*u(j,k); 
				kk++;
			} else if (k == ny-1)
			{
				rhSide(kk,1) = ry*u(j, k-1) + (1.- ry)*u(j, k);
			} else
			{
				rhSide(kk,1) = ry/2.*u(j, k+1) + (1.- ry)*u(j, k) +
							   ry/2.*u(j, k-1);
				kk++;
			}
		}

     return rhSide;
}

// get the matrix B as a packed tridiagonal matrix
Matrix heatADT::StencilB(double &rx)
{
	int i;
	Matrix B(nx, 3);
    for( i=2; i<=(nx-2); i++) {
		B(i,1) = -rx/2.;
		B(i,2) = 1.+ rx;  // Set interior rows
		B(i,3) = -rx/2.;
	}
// Neumann Boundary
	B(1,1)=0.0; B(1,2)= 1.+ rx; B(1,3)=-rx;
	B(nx-1, 1) = -rx; B(nx-1, 2) = 1.+ rx; B(nx-1,3)= 0.0;

	return B;
}

// this function return the right hand side of the second equation
// here u is u^{n+1/2}
Matrix heatADT::getrhsQ(int j,double rx, Matrix &u)
{
	int k,  kk=1;
	Matrix rhSide (ny,1);

	for(k=1; k<ny-1; k++)
	{
		if (j == 1)
		{
			rhSide(kk,1) = rx*u(j+1, k) + (1.- rx)*u(j, k); 
			kk++;
		} else if (j == nx-1)
		{
			rhSide(kk,1) = rx*u(j-1, k) + (1.- rx)*u(j, k);
		} else
		{
			rhSide(kk,1) = rx/2.*u(j+1, k) + (1.- rx)*u(j, k) +
						   rx/2.*u(j-1, k);
			kk++;
		}
	}

     return rhSide;
}

// get the matrix Q2 as a packed tridiagonal matrix
Matrix heatADT::StencilQ(double ry) 
{
	int i=0;
	Matrix B(ny, 3);

    for( i=2; i<=(ny-2); i++) {
		B(i,1) = -ry/2.;
		B(i,2) = 1.+ ry;  // Set interior rows
		B(i,3) = -ry/2.;
	}

// Neumann Boundary
	B(1,1)= 0.; B(1,2)= 1.+ ry; B(1,3)=-ry; 
	B(ny-1, 1) = -ry; B(ny-1, 2) = 1.+ ry; B(ny-1,3)= 0.;

	return B;
}

void heatADT::tridge( Matrix A, Matrix b, Matrix& x) 
{
// Function to solve b = A*x by Gaussian elimination where
// the matrix A is a packed tridiagonal matrix
// Inputs
//   A      Packed tridiagonal matrix, N by N unpacked
//   b      Column vector of length N
// Output 
//   x      Solution of b = A*x; Column vector of length N
// determ   Determinant of A

  // Check that dimensions of a and b are compatible
	int N = A.nRow();
	assert( N == b.nRow() && A.nCol() == 3 );
 
  // Unpack diagonals of triangular matrix into vectors
	Matrix alpha(N,1), beta(N,1), gamma(N,1);
	int i;
	for( i=1; i<=(N-1); i++ ) 
	{
		alpha(i,1) = A(i+1,1);
		beta(i,1) = A(i,2);
		gamma(i,1) = A(i,3);
	}
	beta(N-1,1) = A(N-1,2);

  // Perform forward elimination
	for( i=2; i<=N; i++ )	
	{
		double coeff = alpha(i-1,1)/beta(i-1,1);
		beta(i,1) -= coeff*gamma(i-1,1);
		b(i,1) -= coeff*b(i-1,1);
	}

  // Perform back substitution
	x(N-1,1) = b(N-1,1)/beta(N-1,1);

	for( i=N-1-1; i>=1; i--)
	{ 
		x(i,1) = (b(i,1) - gamma(i,1)*x(i+1,1))/beta(i,1);
	}
}

// solve the equation 
Matrix heatADT::solver(double kapp,
						 double dt, Matrix &u)
{
	int j=1, k=1;
	Matrix x(nx,1), y(ny,1);
	Matrix B(nx, nx), Br(nx,1), Q(ny, ny), Qr(ny,1); 
	Matrix tempUb(nx,ny), tempUq(nx, ny);
	double rx, ry;

	int nn =1;				// divide timestep into nn substep
	double ddt =dt/(double)nn;
	rx = kapp*ddt/(dx*dx); 
	ry = kapp*ddt/(dy*dy);

	for (int kk=0; kk<nn; kk++) 
	{ 
	 for (k=1; k<ny; k++)
	 {
		Br = getrhsB(k, ry, u);
		B = StencilB(rx);
		tridge(B, Br, x); 
	
		for (int i = 1; i<nx; i++)
			tempUb(i,k ) = x(i,1);
	 }

	 for (j=1; j<nx; j++)
	 {
		Qr=getrhsQ(j, rx, tempUb);
		Q=StencilQ(ry);
		tridge(Q, Qr, y);
	
		for (int h=1; h<ny; h++)
			tempUq(j, h) = y(h,1);
	 }
		u = tempUq;
	}

	return tempUq;
}



//  solve Level set equation
Matrix heatADT::FTCSsolver(double mu, double h, double dt, Matrix& u,
						   Image &U0, double &U1, double &U2)
{
    int i, j;
    double r1=1.0, r2=1.0;
	Matrix uNew(nx, ny);
	
	double C1, C2, C3, C4, C, xm, D;

	double epsilon = h/5.; 
	double e = h/100.;

	for(i=2; i<=nx-1; i++)
	for(j=2; j<=ny-1; j++)
	{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-1)),2)/(4.*h*h))+e);
		C2=1./(sqrt(pow((u(i,j) - u(i-1,j)),2)/(h*h) +
				   pow((u(i-1,j+1) - u(i-1,j-1)),2)/(4.*h*h)) +e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-1,j)),2)/(4.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h)) +e);
		C4=1./(sqrt(pow((u(i+1,j-1) - u(i-1,j-1)),2)/(4.*h*h) +
				   pow((u(i,j) - u(i,j-1)),2)/(h*h))+e);
			
		D = 1./PI *epsilon/(epsilon*epsilon +u(i,j)*u(i,j));
		xm = mu*dt *D/(h*h);
		C = 1.+xm*(C1+C2+C3+C4);

		uNew(i,j) = C*u(i,j)+xm*(C1*u(i+1,j) + C2*u(i-1,j) +
									 C3*u(i,j+1) + C4*u(i,j-1)) 
					   + dt*D*(-r1*pow((U0(i,j)-U1),2) + 
					            r2*pow((U0(i,j)-U2),2));

		}

     return uNew;
}


Matrix heatADT::solveDiffu(double kapp, double dt, Matrix& u)
						 //  Image &U0, double &U1, double &U2)
{
    int i, j;
    double h =dx;
	Matrix uNew(nx, ny);
	double xm, V=0.;	
	double C1; //, C2;, C3, C4, C, xm, D;

//	double epsilon = h/5.; 
//	double e = h/100.;
 
	for(i=2; i<=nx-1; i++)
	for(j=2; j<=ny-1; j++)
	{
		xm = kapp*dt/(h*h);
	
		C1 = sqrt(pow((u(i+1,j) - u(i-1,j)),2)/(4.*h*h) +
			       pow((u(i,j+1) - u(i,j-1)),2)/(4.*h*h));
		//if (i == nx/2)
	
	//	cout<<" C1 = "<<C1<<endl;

			V =  C1*20.;//+ pow(C1,2)/9.+ pow(C1,3)/pow(4.,3);
	//	else V= 1.;

		uNew(i,j) = u(i,j)+xm*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) 
							   - 4.*u(i,j)) +  
					   + dt*(V /*u(i,j)*/); 
					            

		}

     return uNew;
}
