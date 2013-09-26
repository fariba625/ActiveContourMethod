/*******************************************************************
 * file:    LevelADT.cpp 
 * version: 1.0
 * project: Project
 * author:  Song Gao
 * date:    July 20, 2003
 * description: Function implementations in class LevelADT
 * 	           
 ******************************************************************/ 

#include <iostream>
#include <fstream>
#include <assert.h>  
#include <math.h>
#include "DMatrix.h"

#include "LevelADT.h"

#define PI 3.1415926

         
// default constructor
LevelADT::LevelADT()
{    
	Image temp;
	nx = temp.getRow(); //256; 
	ny = temp.getCol(); //256; 
	dx = 1.0;
	dy = 1.0;
}

//constructor
LevelADT::LevelADT(int row, int col)
{    
	Image temp(row,col);
	nx = temp.getRow(); //256; 
	ny = temp.getCol(); //256; 
	dx = 1.0;
	dy = 1.0;
}

/*
// Initial conditions 
void LevelADT::init(Matrix &u)
{
	int j,k;

	for(j=1; j<=nx; j++)
	for(k=1; k<=ny; k++)
	{ // a circle with Radius = 80
		u(j, k) = sqrt(pow((j - (double)nx/2.0),2)+ 
				  pow((k - (double)ny/2.0),2))-70.;
	}
}
*/

// this function return the right hand side of the first equation
// here u is u^{n+1/2}
Matrix LevelADT::getrhsB(int k, Matrix &ry1, Matrix &ry2, Matrix &DT, Matrix &u)
{
	int j,  kk=1;
	Matrix rhSide(nx,1);
	
		for(j=1; j<nx-1; j++) 
		{
			if (k == 1)
			{
				rhSide(kk,1) = (ry1(j,k)+ry2(j,k))/2.*u(j, k+1) + 
					           (1.- (ry1(j,k)+ry2(j,k))/2.)*u(j,k) + DT(j,k); 
				kk++;
			} else if (k == ny-1)
			{
				rhSide(kk,1) = (ry1(j,k)+ry2(j,k))/2.*u(j, k-1) + 
							   (1.- (ry1(j,k)+ry2(j,k))/2.)*u(j, k) +DT(j,k);
			} else
			{
				rhSide(kk,1) = ry1(j,k)/2.*u(j, k+1) + 
							   (1.- (ry1(j,k)+ry2(j,k))/2.)*u(j, k) +
							   ry2(j,k)/2.*u(j, k-1) + DT(j,k);
				kk++;
			}
		}

     return rhSide;
}

// get the matrix B as a packed tridiagonal matrix
Matrix LevelADT::StencilB(int j, Matrix &rx1, Matrix &rx2)
{
	int i;
	Matrix B(nx, 3);
    for( i=2; i<=(nx-2); i++) {
		B(i,1) = -rx2(i,j)/2.;
		B(i,2) = 1.+ (rx1(i,j)+rx2(i,j))/2.;  // Set interior rows
		B(i,3) = -rx1(i,j)/2.;
	}
// Neumann Boundary
	B(1,1)=0.0; 
	B(1,2)= 1.+ (rx1(1,j)+rx2(1,j))/2.; 
	B(1,3)= -(rx1(1,j) + rx2(1,j))/2.;
	B(nx-1,1) = -(rx1(nx-1,j) + rx2(nx-1,j))/2.; 
	B(nx-1,2) = 1.+ (rx1(nx-1,j) + rx2(nx-1,j))/2.; 
	B(nx-1,3)= 0.0;

	return B;
}

// this function return the right hand side of the second equation
// here u is u^{n+1/2}
Matrix LevelADT::getrhsQ //(int nx, int ny, int j,double rx, Matrix &u)
	(int j, Matrix &rx1, Matrix &rx2, Matrix &DT, Matrix &u)
{
	int k,  kk=1;
	Matrix rhSide (ny,1);

	for(k=1; k<ny-1; k++)
	{
		if (j == 1)
		{
			rhSide(kk,1) = (rx1(j,k)+rx2(j,k))/2.*u(j+1, k) + 
						   (1.- (rx1(j,k)+rx2(j,k))/2.)*u(j, k) +
						   DT(j,k); 
			kk++;
		} else if (j == nx-1)
		{
			rhSide(kk,1) = (rx1(j,k)+rx2(j,k))/2.*u(j-1, k) + 
						   (1.- (rx1(j,k)+rx2(j,k))/2.)*u(j, k) +
						   DT(j,k);
		} else
		{
			rhSide(kk,1) = rx1(j,k)/2.*u(j+1, k) + 
						   (1.- (rx1(j,k)+rx2(j,k))/2.)*u(j, k) +
						   rx2(j,k)/2.*u(j-1, k) +DT(j,k);
			kk++;
		}
	}

     return rhSide;
}

// get the matrix Q2 as a packed tridiagonal matrix
Matrix LevelADT::StencilQ(int j, Matrix &ry1, Matrix &ry2)
{
	int i;
	Matrix Q(ny, 3);
    for( i=2; i<=(ny-2); i++) {
		Q(i,1) = -ry2(j,i)/2.;
		Q(i,2) = 1.+ (ry1(j,i)+ry2(j,i))/2.;  // Set interior rows
		Q(i,3) = -ry1(j,i)/2.;
	}
// Neumann Boundary
	Q(1,1)=0.0; 
	Q(1,2)= 1.+ (ry1(j,1)+ry2(j,1))/2.; 
	Q(1,3)= -(ry1(1,j) + ry2(j,1))/2.;
	Q(ny-1,1) = -(ry1(j,ny-1) + ry2(j,ny-1))/2.; 
	Q(ny-1,2) = 1.+ (ry1(j,ny-1) + ry2(j,ny-1))/2.; 
	Q(ny-1,3)= 0.0;

	return Q;
}


void LevelADT::tridge( Matrix A, Matrix b, Matrix& x) 
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

// solve the leve equation using ADT 
Matrix LevelADT::solver(double mu, double dt, Matrix &u,Image &U0, 
						double &U1, double &U2)
{
	int i =1, j=1, k=1;
	Matrix x(nx,1), y(ny,1);
	Matrix B(nx, nx), Br(nx,1), Q(ny, ny), Qr(ny,1); 
	Matrix tempUb(nx,ny), tempUq(nx, ny);
	Matrix rx1(nx, ny), rx2(nx, ny), ry1(nx, ny), ry2(nx, ny);
	Matrix DT(nx, ny);

	double r1=1.0, r2=1.0;

	double C1, C2, C3, C4, D;
	double h =dx;
	double epsilon = 20.*dx;// /20.; 
    double e = dx; // /100.;

	for(i=1; i<=nx-1; i++)
	for(j=1; j<=ny-1; j++)
	{
		if (i ==1&&j ==1)
		{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-0)),2)/(1.*h*h)) + e);
		C2=1./(sqrt(pow((u(i,j) - u(i+1,j)),2)/(h*h) +
				   pow((u(i-0,j+1) - u(i-0,j-0)),2)/(1.*h*h)) +e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-0,j)),2)/(1.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h)) +e);
		C4=1./(sqrt(pow((u(i+1,j-0) - u(i-0,j-0)),2)/(1.*h*h) +
				   pow((u(i,j) - u(i,j+1)),2)/(h*h)) +e);
		}else if (i==1&&j> 1) 
		{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-1)),2)/(4.*h*h)) +e);
		C2=1./(sqrt(pow((u(i,j) - u(i+1,j)),2)/(h*h) +
				   pow((u(i-0,j+1) - u(i-0,j-1)),2)/(4.*h*h)) +e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-0,j)),2)/(1.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h)) +e);
		C4=1./(sqrt(pow((u(i+1,j-1) - u(i-0,j-1)),2)/(1.*h*h) +
				   pow((u(i,j) - u(i,j-1)),2)/(h*h))+e);
		} else if (i>1&&j == 1) 
		{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-0)),2)/(1.*h*h)) +e);
		C2=1./(sqrt(pow((u(i,j) - u(i-1,j)),2)/(h*h) +
				   pow((u(i-1,j+1) - u(i-1,j-0)),2)/(1.*h*h)) +e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-1,j)),2)/(4.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h))+e);
		C4=1./(sqrt(pow((u(i+1,j-0) - u(i-1,j-0)),2)/(4.*h*h) +
				   pow((u(i,j) - u(i,j+1)),2)/(h*h))+e);
		}else //(i>1&&j>1)
		{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-1)),2)/(4.*h*h))+e);
		C2=1./(sqrt(pow((u(i,j) - u(i-1,j)),2)/(h*h) +
				   pow((u(i-1,j+1) - u(i-1,j-1)),2)/(4.*h*h))+e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-1,j)),2)/(4.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h))+e);
		C4=1./(sqrt(pow((u(i+1,j-1) - u(i-1,j-1)),2)/(4.*h*h) +
				   pow((u(i,j) - u(i,j-1)),2)/(h*h))+e);
		}
	
		

		D = 1./PI *epsilon/(epsilon*epsilon +u(i,j)*u(i,j));

		DT(i,j) = dt/2.*D*(-r1*pow((U0(i,j)-U1),2) + 
					            r2*pow((U0(i,j)-U2),2) );

		rx1(i,j) = mu*dt/(dx*dx)*D*C1;
		rx2(i,j) = mu*dt/(dx*dx)*D*C2;
		ry1(i,j) = mu*dt/(dy*dy)*D*C3;
		ry2(i,j) = mu*dt/(dy*dy)*D*C4;
	}

	for (k=1; k<ny; k++)
	 {
		Br = getrhsB(k, ry1, ry2, DT, u);
		B = StencilB(k, rx1, rx2);
		tridge(B, Br, x); 
	
		for (int i = 1; i<nx; i++)
			tempUb(i,k ) = x(i,1);
	 }

	 for (j=1; j<nx; j++)
	 {
		Qr=getrhsQ(j, rx1, rx2, DT, tempUb);
		Q=StencilQ(j, ry1, ry2);
		tridge(Q, Qr, y);
	
		for (int h=1; h<ny; h++)
			tempUq(j, h) = y(h,1);
	 }

	return tempUq;
}


//  solve Level set equation
void LevelADT::FTCSsolver(double mu, double dt, Matrix& u,
						   Image &U0, double &U1, double &U2)
{
    int i, j;
	Matrix uNew(nx, ny);

	for(i=1; i<=nx; i++)
	for(j=1; j<=ny; j++)
	{
		if(i==1||i==nx||j==1||j==ny){
			uNew(i,j)=u(i,j);
			continue;
		}
		float ux,uy,uxx,uyy,uxy;
		uy=(u(i,j+1)-u(i,j-1))/2;
		ux=(u(i+1,j)-u(i-1,j))/2;
		float k1,k2,k;
		k2=ux*ux+uy*uy;
		if(k2<1e-3)
			k=0;
		else{
			uxx=(u(i+1,j)-2*u(i,j)+u(i-1,j))/2;
			uyy=(u(i,j+1)-2*u(i,j)+u(i,j-1))/2;
			uxy=(u(i+1,j+1)-u(i+1,j-1)-u(i-1,j+1)+u(i-1,j-1))/4;
			k1=uxx*uy*uy-2*ux*uy*uxy+uyy*ux*ux;
			k=k1/k2;
		}
		float v = u(i,j)+mu*dt*k+dt*(-1.*pow((U0(i,j)-U1),2)+ 
					            pow((U0(i,j)-U2),2));

		uNew(i,j)=v;
	}
	u=uNew;    
}


//  solve Level set equation for phi_1, here phi_1 = u
void LevelADT::phi_1solver(double mu, double dt, Matrix& u,
						   Image &u0, double c11, double c00, 
						   double c10, double c01, Matrix& Hphi2)
{
    int i, j;
	Matrix uNew(nx, ny);	// new phi_1 

	for(i=2; i<=nx-1; i++)
	for(j=2; j<=ny-1; j++)
	{
		if(i==1||i==nx||j==1||j==ny){
			uNew(i,j)=u(i,j);
			continue;
		}
		float ux,uy,uxx,uyy,uxy;
		uy=(u(i,j+1)-u(i,j-1))/2;
		ux=(u(i+1,j)-u(i-1,j))/2;
		float k1,k2,k;
		k2=ux*ux+uy*uy;
		if(k2<1e-3)
			k=0;
		else{
			uxx=(u(i+1,j)-2*u(i,j)+u(i-1,j))/2;
			uyy=(u(i,j+1)-2*u(i,j)+u(i,j-1))/2;
			uxy=(u(i+1,j+1)-u(i+1,j-1)-u(i-1,j+1)+u(i-1,j-1))/4;
			k1=uxx*uy*uy-2*ux*uy*uxy+uyy*ux*ux;
			k=k1/k2;
		}
		uNew(i,j) =u(i,j)+mu*dt*k
					   + dt*(-pow((u0(i,j)-c11),2)*Hphi2(i,j) 
							   -pow((u0(i,j)-c10),2)*(1-Hphi2(i,j))
							   +pow((u0(i,j)-c01),2)*Hphi2(i,j) //+
					           +pow((u0(i,j)-c00),2)*(1-Hphi2(i,j)));
	}
	u=uNew;	
}


//  solve Level set equation for phi_2, here phi_2 = u
void LevelADT::phi_2solver(double mu, double dt, Matrix& u,
						   Image &u0, double c11, double c00, 
						   double c10, double c01, Matrix& Hphi1)
{
    int i, j;
	Matrix uNew(nx, ny);	// new phi_2[i][j] 
	
	for(i=2; i<=nx-1; i++)
	for(j=2; j<=ny-1; j++)
	{
		if(i==1||i==nx||j==1||j==ny){
			uNew(i,j)=u(i,j);
			continue;
		}
		float ux,uy,uxx,uyy,uxy;
		uy=(u(i,j+1)-u(i,j-1))/2;
		ux=(u(i+1,j)-u(i-1,j))/2;
		float k1,k2,k;
		k2=ux*ux+uy*uy;
		if(k2<1e-3)
			k=0;
		else{
			uxx=(u(i+1,j)-2*u(i,j)+u(i-1,j))/2;
			uyy=(u(i,j+1)-2*u(i,j)+u(i,j-1))/2;
			uxy=(u(i+1,j+1)-u(i+1,j-1)-u(i-1,j+1)+u(i-1,j-1))/4;
			k1=uxx*uy*uy-2*ux*uy*uxy+uyy*ux*ux;
			k=k1/k2;
		}
		uNew(i,j) = u(i,j)+mu*dt*k
					   + dt*(-pow((u0(i,j)-c11),2)*Hphi1(i,j) 
							   +pow((u0(i,j)-c10),2)*(Hphi1(i,j))
							   -pow((u0(i,j)-c01),2)*(1-Hphi1(i,j))
					           +pow((u0(i,j)-c00),2)*(1-Hphi1(i,j)));
		}
     u=uNew;
}


// solve the leve equation using ADT For phi_1 
Matrix LevelADT::ADTsolver1(double mu, double dt, Matrix &u,Image &u0, 
							double c11, double c00, 
						   double c10, double c01, Matrix& Hphi2)
{
	int i =1, j=1, k=1;
	Matrix x(nx,1), y(ny,1);
	Matrix B(nx, nx), Br(nx,1), Q(ny, ny), Qr(ny,1); 
	Matrix tempUb(nx,ny), tempUq(nx, ny);
	Matrix rx1(nx, ny), rx2(nx, ny), ry1(nx, ny), ry2(nx, ny);
	Matrix DT(nx, ny);

	//double r1=1.0, r2=1.0;

	double C1, C2, C3, C4, D;
	double h =dx;
	double epsilon = 15.*dx;// /20.; 
    double e = dx/100.;

	for(i=1; i<=nx-1; i++)
	for(j=1; j<=ny-1; j++)
	{
		if (i ==1&&j ==1)
		{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-0)),2)/(1.*h*h)) + e);
		C2=1./(sqrt(pow((u(i,j) - u(i+1,j)),2)/(h*h) +
				   pow((u(i-0,j+1) - u(i-0,j-0)),2)/(1.*h*h)) +e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-0,j)),2)/(1.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h)) +e);
		C4=1./(sqrt(pow((u(i+1,j-0) - u(i-0,j-0)),2)/(1.*h*h) +
				   pow((u(i,j) - u(i,j+1)),2)/(h*h)) +e);
		}else if (i==1&&j> 1) 
		{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-1)),2)/(4.*h*h)) +e);
		C2=1./(sqrt(pow((u(i,j) - u(i+1,j)),2)/(h*h) +
				   pow((u(i-0,j+1) - u(i-0,j-1)),2)/(4.*h*h)) +e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-0,j)),2)/(1.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h)) +e);
		C4=1./(sqrt(pow((u(i+1,j-1) - u(i-0,j-1)),2)/(1.*h*h) +
				   pow((u(i,j) - u(i,j-1)),2)/(h*h))+e);
		} else if (i>1&&j == 1) 
		{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-0)),2)/(1.*h*h)) +e);
		C2=1./(sqrt(pow((u(i,j) - u(i-1,j)),2)/(h*h) +
				   pow((u(i-1,j+1) - u(i-1,j-0)),2)/(1.*h*h)) +e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-1,j)),2)/(4.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h))+e);
		C4=1./(sqrt(pow((u(i+1,j-0) - u(i-1,j-0)),2)/(4.*h*h) +
				   pow((u(i,j) - u(i,j+1)),2)/(h*h))+e);
		}else //(i>1&&j>1)
		{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-1)),2)/(4.*h*h))+e);
		C2=1./(sqrt(pow((u(i,j) - u(i-1,j)),2)/(h*h) +
				   pow((u(i-1,j+1) - u(i-1,j-1)),2)/(4.*h*h))+e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-1,j)),2)/(4.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h))+e);
		C4=1./(sqrt(pow((u(i+1,j-1) - u(i-1,j-1)),2)/(4.*h*h) +
				   pow((u(i,j) - u(i,j-1)),2)/(h*h))+e);
		}
	
		

		D = 1./PI *epsilon/(epsilon*epsilon +u(i,j)*u(i,j));

		DT(i,j) = dt/2.*D*(-pow((u0(i,j)-c11),2)*Hphi2(i,j) 
							   -pow((u0(i,j)-c10),2)*(1-Hphi2(i,j))
							   +pow((u0(i,j)-c01),2)*Hphi2(i,j)
					           +pow((u0(i,j)-c00),2)*(1-Hphi2(i,j)));

		rx1(i,j) = mu*dt/(dx*dx)*D*C1;
		rx2(i,j) = mu*dt/(dx*dx)*D*C2;
		ry1(i,j) = mu*dt/(dy*dy)*D*C3;
		ry2(i,j) = mu*dt/(dy*dy)*D*C4;
	}

	for (k=1; k<ny; k++)
	 {
		Br = getrhsB(k, ry1, ry2, DT, u);
		B = StencilB(k, rx1, rx2);
		tridge(B, Br, x); 
	
		for (int i = 1; i<nx; i++)
			tempUb(i,k ) = x(i,1);
	 }

	 for (j=1; j<nx; j++)
	 {
		Qr=getrhsQ(j, rx1, rx2, DT, tempUb);
		Q=StencilQ(j, ry1, ry2);
		tridge(Q, Qr, y);
	
		for (int h=1; h<ny; h++)
			tempUq(j, h) = y(h,1);
	 }

	return tempUq;
}


// solve the leve equation using ADT For phi_1 
Matrix LevelADT::ADTsolver2(double mu, double dt, Matrix &u,Image &u0, 
							double c11, double c00, 
						   double c10, double c01, Matrix& Hphi1)
{
	int i =1, j=1, k=1;
	Matrix x(nx,1), y(ny,1);
	Matrix B(nx, nx), Br(nx,1), Q(ny, ny), Qr(ny,1); 
	Matrix tempUb(nx,ny), tempUq(nx, ny);
	Matrix rx1(nx, ny), rx2(nx, ny), ry1(nx, ny), ry2(nx, ny);
	Matrix DT(nx, ny);

	//double r1=1.0, r2=1.0;

	double C1, C2, C3, C4, D;
	double h =dx;
	double epsilon = 15.*dx;// /20.; 
    double e = dx/100.;

	for(i=1; i<=nx-1; i++)
	for(j=1; j<=ny-1; j++)
	{
		if (i ==1&&j ==1)
		{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-0)),2)/(1.*h*h)) + e);
		C2=1./(sqrt(pow((u(i,j) - u(i+1,j)),2)/(h*h) +
				   pow((u(i-0,j+1) - u(i-0,j-0)),2)/(1.*h*h)) +e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-0,j)),2)/(1.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h)) +e);
		C4=1./(sqrt(pow((u(i+1,j-0) - u(i-0,j-0)),2)/(1.*h*h) +
				   pow((u(i,j) - u(i,j+1)),2)/(h*h)) +e);
		}else if (i==1&&j> 1) 
		{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-1)),2)/(4.*h*h)) +e);
		C2=1./(sqrt(pow((u(i,j) - u(i+1,j)),2)/(h*h) +
				   pow((u(i-0,j+1) - u(i-0,j-1)),2)/(4.*h*h)) +e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-0,j)),2)/(1.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h)) +e);
		C4=1./(sqrt(pow((u(i+1,j-1) - u(i-0,j-1)),2)/(1.*h*h) +
				   pow((u(i,j) - u(i,j-1)),2)/(h*h))+e);
		} else if (i>1&&j == 1) 
		{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-0)),2)/(1.*h*h)) +e);
		C2=1./(sqrt(pow((u(i,j) - u(i-1,j)),2)/(h*h) +
				   pow((u(i-1,j+1) - u(i-1,j-0)),2)/(1.*h*h)) +e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-1,j)),2)/(4.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h))+e);
		C4=1./(sqrt(pow((u(i+1,j-0) - u(i-1,j-0)),2)/(4.*h*h) +
				   pow((u(i,j) - u(i,j+1)),2)/(h*h))+e);
		}else //(i>1&&j>1)
		{
		C1=1./(sqrt(pow((u(i+1,j) - u(i,j)),2)/(h*h) +
			       pow((u(i,j+1) - u(i,j-1)),2)/(4.*h*h))+e);
		C2=1./(sqrt(pow((u(i,j) - u(i-1,j)),2)/(h*h) +
				   pow((u(i-1,j+1) - u(i-1,j-1)),2)/(4.*h*h))+e);
		C3=1./(sqrt(pow((u(i+1,j) - u(i-1,j)),2)/(4.*h*h) +
				   pow((u(i,j+1) - u(i,j)),2)/(h*h))+e);
		C4=1./(sqrt(pow((u(i+1,j-1) - u(i-1,j-1)),2)/(4.*h*h) +
				   pow((u(i,j) - u(i,j-1)),2)/(h*h))+e);
		}
	
		

		D = 1./PI *epsilon/(epsilon*epsilon +u(i,j)*u(i,j));

		DT(i,j) = dt/2.*D*(-pow((u0(i,j)-c11),2)*Hphi1(i,j) 
							   +pow((u0(i,j)-c10),2)*(Hphi1(i,j))
							   -pow((u0(i,j)-c01),2)*(1-Hphi1(i,j))
					           +pow((u0(i,j)-c00),2)*(1-Hphi1(i,j)));

		rx1(i,j) = mu*dt/(dx*dx)*D*C1;
		rx2(i,j) = mu*dt/(dx*dx)*D*C2;
		ry1(i,j) = mu*dt/(dy*dy)*D*C3;
		ry2(i,j) = mu*dt/(dy*dy)*D*C4;
	}

	for (k=1; k<ny; k++)
	 {
		Br = getrhsB(k, ry1, ry2, DT, u);
		B = StencilB(k, rx1, rx2);
		tridge(B, Br, x); 
	
		for (int i = 1; i<nx; i++)
			tempUb(i,k ) = x(i,1);
	 }

	 for (j=1; j<nx; j++)
	 {
		Qr=getrhsQ(j, rx1, rx2, DT, tempUb);
		Q=StencilQ(j, ry1, ry2);
		tridge(Q, Qr, y);
	
		for (int h=1; h<ny; h++)
			tempUq(j, h) = y(h,1);
	 }

	return tempUq;
}

