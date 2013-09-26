
#ifndef MATRIX_CLASS
#define MATRIX_CLASS

#include <assert.h>  // Defines the assert function.

class Matrix
{
public:
	Matrix ();	// Default Constructor
 	Matrix(int, int);	// Regular Constructor
 	Matrix(const Matrix& mat);	// Copy Constructor.
	~Matrix();	// Destructor
	Matrix& operator=(const Matrix& mat);	// Assignment operator
	int nRow() const;	// Simple "get" functions.
	int nCol() const;
 	double& operator() (int, int);	// Parenthesis operator
 // Parenthesis operator function (const version).
	const double& operator() (int, int) const;
// Set function. Sets all elements of a matrix to a given value.
	void set(double value);
	void writeFile(char *fname);

private:

// Matrix data.
	int nRow_, nCol_;  // Number of rows, columns
	double* data_;     // Pointer used to allocate memory for data.

// Private copy function.
// Copies values from one Matrix object to another.
	void copy(const Matrix& mat);
};

#endif	// Class Matrix
