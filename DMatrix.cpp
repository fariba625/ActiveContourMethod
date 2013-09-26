
#include "DMatrix.h"  // Defines the assert function.
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Default Constructor. Creates a 1 by 1 matrix; sets value to zero.
Matrix::Matrix ()
{
	nRow_ = 256; nCol_ = 256;
	data_ = new double [256*256];  // Allocate memory
	set(0.0);                // Set value of data_[0] to 0.0
}

// Regular Constructor. Creates an nR by nC matrix; sets values to zero.
// If number of columns is not specified, it is set to 1.
Matrix::Matrix(int nR, int nC = 1)
{
	assert(nR > 0 && nC > 0);    // Check that nC and nR both > 0.
	nRow_ = nR; nCol_ = nC;
	data_ = new double [nR*nC];  // Allocate memory
	assert(data_ != 0);          // Check that memory was allocated
	set(0.0);                    // Set values of data_[] to 0.0
}

// Copy Constructor.
// Used when a copy of an object is produced
// (e.g., passing to a function by value)
Matrix::Matrix(const Matrix& mat)
{
	this->copy(mat);   // Call private copy function.
}

// Destructor. Called when a Matrix object goes out of scope.
Matrix::~Matrix()
{
	delete [] data_;   // Release allocated memory
}

// Assignment operator function.
// Overloads the equal sign operator to work with
// Matrix objects.
Matrix& Matrix::operator=(const Matrix& mat)
{
	if( this == &mat ) return *this;  // If two sides equal, do nothing.
	delete [] data_;                  // Delete data on left hand side
	this->copy(mat);                  // Copy right hand side to l.h.s.
	return *this;
}

// Simple "get" functions. Return number of rows or columns.
int Matrix::nRow() const
{
	return nRow_;
}

int Matrix::nCol() const
{
	return nCol_;
}

// Parenthesis operator function.
// Allows access to values of Matrix via (i,j) pair.
// Example: a(1,1) = 2*b(2,3);
// If column is unspecified, take as 1.
double& Matrix::operator() (int i, int j = 1)
{
	assert(i > 0 && i <= nRow_);          // Bounds checking for rows
	assert(j > 0 && j <= nCol_);          // Bounds checking for columns
	return data_[ nCol_*(i-1) + (j-1) ];  // Access appropriate value
}

// Parenthesis operator function (const version).
const double& Matrix::operator() (int i, int j = 1) const
{
	assert(i >= 0 && i < nRow_);          // Bounds checking for rows
	assert(j >= 0 && j < nCol_);          // Bounds checking for columns
	return data_[ nCol_*(i-1) + (j-1) ];  // Access appropriate value
}

// Set function. Sets all elements of a matrix to a given value.
void Matrix::set(double value)
{
	int i, iData = nRow_*nCol_;
	for( i=0; i<iData; i++ )
		data_[i] = value;
}

void Matrix::copy(const Matrix& mat)
{
	nRow_ = mat.nRow_;
	nCol_ = mat.nCol_;
	int i, iData = nRow_*nCol_;
	data_ = new double [iData];
	for(i = 0; i<iData; i++ )
		data_[i] = mat.data_[i];
}

void Matrix::writeFile(char *fname)
{
	ofstream ofp;
	int i, j;

	ofp.open(fname, ios::out);

	if (!ofp) {
		cout << "Can't write image: " << fname << endl;
		exit(1);
	}

	for (i=1; i<nRow_-1; i++)
	{
	  for (j=1; j<nCol_-1; j++)
		ofp<<data_[ nCol_*(i-1) + (j-1) ] <<endl;
	//	ofp<<endl;
	}

  ofp.close();
}
