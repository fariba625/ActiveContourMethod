
#include "cv.h"
#include "highgui.h"

#include "LevelADT.h"
#include "2DheatADT.h"

#ifndef Image_CLASS
#define Image_CLASS

#define GRAY     10
#define BINARY   11

class Image
{
 public:
  // constructors and destructor
	Image();                          // default constructor
	Image(int, int);                   // constructor with row, column

	Image(Image &);                  // copy constructor
	~Image();                         // destructor

  // create a Image
	void createImage();
	void createImage(int,             // row
		    int);            // column

  // get and set functions
	int getRow() const;                // get row number / the number of sample
	int getCol() const;                // get column number / the number of feature
	int getMaximum() const;              // get the maximum pixel value
	int getMinimum() const;              // get the mininum pixel value
	void setRow(int);                  // set row number
	void setCol(int);                  // set column number

  // read in data to a Image
	void readImage(char *);
	void readPGMRAW(char *);           // read the image into a matrix (row,col)
	void writePGMRAW(char *);

  // operator overloading functions
	double & operator()(int,           // row index
		    int);            // column index
	Image operator=(Image);          // = operator overloading
	Image operator+(Image);          // overloading + operator
	Image operator-(Image);          // overloading - operator
	Image operator*(Image);          // overloading * operator (element-wised multiplication)
	Image operator->*(Image);        // overloading ->* operator (matrix multiplication)
	Image operator/(double&);            // overloading scale division

  // Image manipulations
	Image subImage(int,              // partial Image (start row index)
		   int,              // start column index
		   int,              // end row index
		   int);             // end column index

	void allocateMem();
	void deAllocateMem();
	Image threshold(double, int nt = BINARY); // thresholding an image into gray or binary
	Image Normalization(int, int);
	Image negative();  //  compute the image negative
	double getPixel(int);

//	Image InitialCurve();

//	Image moveCurve();

//	void initializePhi();
	Image creatImage(double R1, double R2);
	Image initializePhiPart_2(double, double);
	Image initializePhiPart_4(double, double);
	Image initializePhiPart_6(double, double);
	Image initializePhiPart_m(double, double);
	Image initializePhiPart_n(double, double);
	Image initializePhiPart_s(int, int);

	Image updateD();
	Image DeNoise(Image&);

	Image updateLevel(double&);
	void updateLevel1(double&);
	double energy1();
 private:
	double *pixel;                     // Image buffer
	int row;                           // number of rows
	int col;                           // number of columns
	double **phi1;						// level set function
	double **phi2;						// level set function

	void outofMemory();
	void average(double &c1, double &c0);
	void average(double &c11, double &c00,double &c10, double &c01);
	void average1(double &c11, double &c00,double &c10, double &c01);
	double energyDifference(int x, int y, double&, double&, double&, double&);
	double energyDifference2(int x, int y, double&, double&, double&, double&);
	void average(double&, double&, double&, double&, double&, double&, double&, double&);

	void trans(IplImage*iout);
	void trans1(IplImage*iout);
	void smooth();

//	double minDistance(Matrix &A, Matrix &B);
//	double normal(Matrix &A);
};


#endif










