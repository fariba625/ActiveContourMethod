


#include "Image.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define a1 1

#define PI 3.1415926
#define ROTATE(a,i,j,k,l) g=a(i,j);h=a(k,l);a(i,j)=g-s*(h+g*tau);\
   a(k,l)=h+s*(g-h*tau);
#define ROTATEM(a,i,j,k,l) g=a[i*col+j];h=a[k*col+l];a[i*col+j]=g-s*(h+g*tau);\
   a[k*col+l]=h+s*(g-h*tau);

extern float mu;
static int compareIntegers (const void *e1, const void *e2)
{
	return (*((int *) e1) - *((int *) e2));
}


double min3(double V1, double V2, double V3)
{
		return ((V1<=V2)?((V1<=V3)?V1:V3):(V2<=V3)?V2:V3);
}

// default constructor
Image::Image()
{
	createImage(256,256);//(324, 256);
}

// constructor when knowing row and column
Image::Image(int nr, int nc=1)
{
	createImage(nr, nc);
}

// copy constructor
Image::Image(Image &m)
{
	int i, j;

	createImage(m.getRow(), m.getCol());      // allocate memory

	for (i=0; i<row; i++)
    for (j=0; j<col; j++)
		pixel[i * col + j] = m(i, j);
}

void Image::createImage()
{
	int i;

	pixel = (double *) new double [(row+1) * (col+1)];
	if (!pixel)
	outofMemory();

	for (i=0; i<(row+1)*(col+1); i++)
		pixel[i] = 0.0;
}

// allocate memory for the image
void Image::createImage(int nr, int nc)
{
	int i;

	row = nr;
	col = nc;

	pixel = (double *) new double [(row+1) * (col+1)];

	if (!pixel)
		outofMemory();

	for (i=0; i<(row+1)*(col+1); i++)
		pixel[i] = 0.0;

	allocateMem ( );
}

// output out of memory error
void Image::outofMemory()
{
	cerr << "Out of memory!\n";
	exit(1);
}

// destructor
Image::~Image()
{
	if (pixel)
    delete [] pixel;       // free the pixel buffer

	deAllocateMem ( );

}

// get number of rows
int Image::getRow() const
{
	return row;
}

// get number of columns
int Image::getCol() const
{
	return col;
}

// set number of rows
void Image::setRow(int nr)
{
	row = nr;
}

// set number of columns
void Image::setCol(int nc)
{
	col = nc;
}

// get the maximum pixel value of the image
int Image::getMaximum() const
{
	int tmax = 0;
	int i, j; // k;

	for (i=0; i<row; i++)
    for (j=0; j<col; j++)
		if (tmax < pixel[i*col+j])
			tmax = (int)pixel[i*col+j];

  return tmax;
}

// get the minimum pixel value of the image
int Image::getMinimum() const
{
	double tmin = 255.0;
	int i, j;

	for (i=0; i<row; i++)
    for (j=0; j<col; j++)
		if (tmin  > pixel[i*col+j])
			tmin = pixel[i*col+j];

  return (int)tmin;
}


// read in data to a matrix
void Image::readImage(char *fname)
{
	ifstream infile;
	int i, j;

	infile.open(fname, ios::in);

	if (!infile)
	{
		cerr << "Can't read data file: " << fname << endl;
		exit(1);
	}

	createImage(row, col);

	for (i=0; i<row; i++)
    for (j=0; j<col; j++)
		infile >> pixel[i*col+j];
}

// overloading () operator
double & Image::operator()(int i, int j=1)
{
  return pixel[i * col + j];
}

// overloading = operator
Image Image::operator=(Image m)
{
	int i, j;

	createImage(m.getRow(), m.getCol());             // allocate memory

	for (i=0; i<row; i++)
    for (j=0; j<col; j++)
		pixel[i * col + j] = m(i, j);

	return *this;
}

// overloading + operator
Image Image::operator+(Image m)
{
	int i, j, nr, nc;
	Image temp(row, col);

	nr = m.getRow();
	nc = m.getCol();

	if (nr != row || nc != col)
	{
		cerr << "Matrices are not of the same size, cannot do addition\n";
		exit(3);
	}

//	temp.createImage(row, col);

	for (i=0; i<row; i++)
    for (j=0; j<col; j++)
	{
		if ((pixel[i*col+j] + m(i, j))>= 255)
			temp(i, j) = 255;
		else if ((pixel[i*col+j] + m(i, j)) <=0)
			temp(i, j) = 0;
		else temp(i, j)=pixel[i*col+j] + m(i, j);
	}

	return temp;
}

// overloading - operator
Image Image::operator-(Image m)
{
	int i, j, nr, nc;
	Image temp;

	nr = m.getRow();
	nc = m.getCol();

	if (nr != row || nc != col)
	{
		cerr << "Matrices are not of the same size, cannot do subtraction\n";
		exit(3);
	}

	temp.createImage(row, col);

	for (i=0; i<row; i++)
    for (j=0; j<col; j++)
		temp(i, j) = pixel[i*col+j] - m(i, j);

	return temp;
}

// overloading * operator (piecewise multiplication)
Image Image::operator*(Image m)
{
	int i, j, nr, nc;
	Image temp;

	nr = m.getRow();
	nc = m.getCol();

	if (nr != row || nc != col)
	{
		cerr << "Matrices are not of the same size, cannot do multiplication\n";
		exit(3);
	}

	temp.createImage(row, col);

	for (i=0; i<row; i++)
    for (j=0; j<col; j++)
		temp(i, j) = pixel[i*col+j] * m(i, j);

	return temp;
}

// overloading / operator
Image Image::operator/(double &factor)
{
	int i, j;
	Image temp;

	if (factor == 0.0)
	{
		cout << "Image diveded by 0, can't do division\n";
		exit(1);
	}

	temp.createImage(row, col);

	for (i=0; i<row; i++)
    for (j=0; j<col; j++)
        temp(i, j) = pixel[(i*col+j)]/factor;

	return temp;
}


// overloading ->* operator (matrix multiplication)
Image Image::operator->*(Image m)
{
	int i, j, nr, nc, p;
	Image temp;
	double tmp;

	nr = m.getRow();
	nc = m.getCol();

	if (col != nr)
	{
		cerr << "Image size is not consistent, cannot do multiplication\n";
		exit(3);
	}

	temp.createImage(row, nc);

	for (i=0; i<row; i++)
    for (j=0; j<nc; j++)
	{
      tmp = 0.0;
      for (p=0; p<col; p++)
	  {
		tmp += pixel[i*col+p] * m(p, j);
		temp(i, j) = tmp;
      }
    }

  return temp;
}


// return partial matrix
Image Image::subImage(int startrow, int startcol, int endrow, int endcol)
{
	int i, j;
	Image temp;

	if (startrow > endrow || startcol > endcol || endrow > row -1 || endcol > col - 1)
	{
		cerr << "Check size\n";
		exit(1);
	}

	temp.createImage(endrow-startrow+1, endcol-startcol+1);

	for (i=startrow; i<=endrow; i++)
    for (j=startcol; j<=endcol; j++)
      temp(i-startrow, j-startcol) = pixel[i*col+j];

	return temp;
}


// read PGM image from a file
void Image::readPGMRAW(char *fname)
{
	FILE *pgm;
	char line[512], intbuf[100], ch;
	int type, nc, nr, maxval, i, j, k, found;
    unsigned char *img;


	// *.PGM File Open
	if ((pgm = fopen(fname, "rb")) == NULL)
		cout << "\n		Image " << fname << " can not be opened correctly!" << endl;

	//  Scan pnm type information, expecting P5
	fgets(line, 511, pgm);
	if (line[0]=='#') fgets(line, 511, pgm);
	sscanf(line, "P%d", &type);
	if (type != 5 && type != 2) {
		fclose(pgm);
	}

	//Get dimensions of pgm
	fgets(line, 511, pgm);
	if (line[0]=='#') fgets(line, 511, pgm);
	sscanf(line, "%d %d", &nc, &nr);

	// Get maxval
	fgets(line, 511, pgm);
	if (line[0]=='#') fgets(line, 511, pgm);
	sscanf(line, "%d", &maxval);
	if (maxval > 255) {
		fclose(pgm);
	}

	row = nr;
	col = nc;
	// Memory allocation for input_image
	img = new unsigned char[nr*nc] ;

	if (type == 5)
		fread(img, 1, nr*nc, pgm);
	else if (type == 2) {
		for (i = 0; i < nr; i++) {
			for (j = 0; j < nc; j++) {
				k = 0;  found = 0;
				while (!found) {
					ch = (char) fgetc(pgm);
					if (ch >= '0' && ch <= '9') {
						intbuf[k] = ch;  k++;
					} else {
						if (k != 0) {
							intbuf[k] = '\0';
							found = 1;
						}
					}
				}
				img[i*nc+j] = atoi(intbuf);
			}
		}
	} else {
		fclose(pgm);
	}

	// File Close
	fclose(pgm);

	createImage(row, col);

	for (i=0; i<row; i++)
		for (j=0; j<col; j++)
			pixel[i*col+j] = (double)img[i*col+j];
}



// write to an image
void Image::writePGMRAW(char *fname)
{
	ofstream ofp;
	int i, j;
	unsigned char *img;

	ofp.open(fname, ios::out);

	if (!ofp)
		cerr << "Can't write image: " << fname << endl;

  // Write the format ID
	ofp << "P5" << endl;
	ofp << col << " " << row << endl;
	ofp << 255 << endl;

  // convert the image data type back to unsigned char
	img = (unsigned char *) new unsigned char [row * col];

	for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      img[i*col+j] = (unsigned char)pixel[i*col+j];

	ofp.write(img, (row * col * sizeof(unsigned char)));

	ofp.close();
	delete []img;
}


// compute the negative of the image
Image Image::negative()
{
  Image temp;
  int i, j;

  temp.setRow(row);
  temp.setCol(col);
//  temp.setType(type);
  temp.createImage();

  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
	temp(i, j) = 255 - pixel[i*col + j];

  return temp;
}


// thresholding image to a gray-level or a binary image
Image Image::threshold(double thresh, int nt)
{
	Image temp;
	int i, j;

	temp.setRow(row);
	temp.setCol(col);
	temp.createImage();

	if (nt != GRAY && nt != BINARY)
		nt = BINARY;

	for (i=0; i<row; i++)
	for (j=0; j<col; j++)
	{
		if (pixel[i*col+j] >= thresh)
		   temp(i, j) = (nt == GRAY) ? pixel[i*col+j] : 255;
		else
		temp(i, j) = 0;//pixel[i*col+j];//0;
	}

  return temp;
}


// linear normalization
Image Image::Normalization(int MM, int NN)
{
	int i=0, j=0, m=0, n=0;
	Image temp;

	temp.setRow(MM);
	temp.setCol(NN);
	temp.createImage();

	for (i=0; i<row; i++)
		for (j=0; j<col; j++) {
			if(pixel[i*col+j]>0) {
				m = (int)floor(i*MM/row);
				n = (int)floor(j*NN/col);
				temp(m,n)= pixel[i*col+j];
			}
		  }

			return temp;

}


double Image::getPixel(int i)
{
	return pixel[i];
}

// initial a curve in an image according to the initial condition of the PDE
/*
Image Image::InitialCurve()
{
	int i, j;
	double e =1.;

	Image temp (row, col);

	LevelADT tempX;
	Matrix phi(row, col);

	tempX.init(phi);


	for (i=0; i<row; i++)
	for (j=0; j<col; j++)
	{
		if (phi(i+1, j+1) >= -e&& phi(i+1, j+1)<= e)
			temp(i,j)= 155;
		else
			temp(i,j)= 0;
	}

	return temp;
}


Image Image::moveCurve()
{
	int i, j;
	double e =1./8.;

	Image temp (row, col);

	LevelADT X1;
	Matrix phi(row, col);

	int nx, ny;
	double dx, dy, kapp, dt;

	nx = 256;
	ny = 256;
	dx = 1.;
	dy = 1.;
	kapp =14.0;
	dt = 1.5;
	Matrix uR(nx, ny);

	X1.init(uR);

	for (int tstep = 0; tstep<50; tstep++)
	{
	;//	uR	= X1.solver(nx, ny, dx, dy, kapp, dt, uR);
	}
	phi = uR;

	for (i=0; i<row; i++)
	for (j=0; j<col; j++)
	{
		if (phi(i+1, j+1) >-e && phi(i+1, j+1)<= e)
		{
			temp(i,j)= 155;
		}else
			temp(i,j)= 0;
	}

	return temp;
}
*/

// Allocate memory
void Image::allocateMem( )
{
	phi1 = new double *[row];

	for (int i = 0; i < row; i++)
	{
 	 phi1[i] = new double [col];
 	 if (!phi1[i])
		{
	 	 cerr << "Image: Not enough memory for allocateMem!!"
			<< endl;
	 	 exit(-1);
		}
	 for (int j = 0; j < col; j++)
		phi1[i][j] = 10;
	}
	// phi2
		phi2 = new double *[row];

	for (i = 0; i < row; i++)
	{
 	 phi2[i] = new double [col];
 	 if (!phi2[i])
		{
	 	 cerr << "Image: Not enough memory for allocateMem!!"
			<< endl;
	 	 exit(-1);
		}
	 for (int j = 0; j < col; j++)
		phi2[i][j] = 10;
	}
}

void Image::deAllocateMem ( )
{
	 for (int i = 0; i < row; i++)
	 {
		 delete phi1[i];
		 delete phi2[i];
	 }
	 delete phi1;
	 delete phi2;
}

//initial a curve acrroding phi values (square)
Image Image::initializePhiPart_s(int L, int W)
{
	int i, j;
	double x0 = (double)row;	//cneter of the circle
	double y0 = (double)col;
	double dx =50.;
	double dis1, dis4;
	double e =1.;

//	initializePhi();
	Image temp(row, col);
	allocateMem();

	for (i=0; i<row; i++)
	for (j=0; j<col; j++)
	{
		temp(i,j)= 0;
		if(abs(i-row/2)<L/2&&abs(j-col/2)<W/2){
//		if(sqrt((i-row/2)*(i-row/2)+(j-col/2)*(j-col/2))<L/2){
			phi1[i][j]= 1;
			phi2[i][j]= 1;
		}
		else{
			phi1[i][j]= -1;
			phi2[i][j]= -1;
		}


	}

		return temp;
}



//initial a curve acrroding phi values (2 circles)
Image Image::initializePhiPart_2(double R1, double R2)
{
	int i, j;
	double x0 = (double)row;	//cneter of the circle
	double y0 = (double)col;
	double dx =50.;
	double dis1, dis4;
	double e =1.;

//	initializePhi();
	Image temp(row, col);
	allocateMem();

	for (i=0; i<row; i++)
	for (j=0; j<col; j++)
	{
		dis1 = pow((i - (x0/3.+ dx)),2)+ pow((j - 1.5*y0/3.),2);
		dis4 = pow((i - (2.*x0/3.-dx)),2)+ pow((j - 1.5*y0/3.),2);

		if(sqrt(dis1) > R1-e && sqrt(dis1) <= R1+e )
			temp(i,j)= 255;
		else if (sqrt(dis4) > R2-e && sqrt(dis4) <= R2+e )
			temp(i,j)= -255;
		else
			temp(i,j)= 0;

		if(sqrt(dis1)<=R1)
			phi1[i][j]= 1;
		else if (sqrt(dis1) > R1)
			phi1[i][j]= -1;
//		else phi1[i][j]= 0;

		if(sqrt(dis4)<=R2-0)
			phi2[i][j]= 1;
		else if (sqrt(dis4) > R2+0)
			phi2[i][j]= -1;
//		phi2[i][j]= 1;
	}

		return temp;
}

//initial a curve acrroding phi values (4 circles)
Image Image::initializePhiPart_4(double R1, double R2)
{
	int i, j;
	double x0 = (double)row;	//cneter of the circle
	double y0 = (double)col;
	double dx = 30.;
	double dis1, dis2, dis3, dis4;
	double e =0.5;

//	initializePhi();
	Image temp(row, col);
	allocateMem();

	for (i=0; i<row; i++)
	for (j=0; j<col; j++)
	{
		dis1 = pow((i - (x0/4.+dx)),2)+ pow((j - y0/4.-10),2);
		dis3 = pow((i - (x0/4.+dx)),2)+ pow((j - 3.*y0/4.+10),2);
		dis2 = pow((i - (3.*x0/4.- dx+20)),2)+ pow((j - y0/4.-10.),2);
		dis4 = pow((i - (3.*x0/4.- dx+20)),2)+ pow((j - 3.*y0/4.+ 10.),2);

		if(sqrt(dis1) > R1-e && sqrt(dis1) <= R1+e )
			temp(i,j)= 255;
		else if (sqrt(dis2) > R2-e && sqrt(dis2) <= R2+e )
			temp(i,j)= -255;
		else if (sqrt(dis3) > R1-e && sqrt(dis3) <= R1+e )
			temp(i,j)= 255;
		else if (sqrt(dis4) > R2-e && sqrt(dis4) <= R2+e )
			temp(i,j)= -255;
		else
			temp(i,j)= 0;

		if(sqrt(dis1)<R1-e||sqrt(dis3)<R1-e)
			phi1[i][j]= 1;
		else
			phi1[i][j]= -1;

		if(sqrt(dis2)<R2-0||sqrt(dis4)<R2-0)
			phi2[i][j]= 1;
		else
			phi2[i][j]= -1;
	}

		return temp;
}

//initial a curve acrroding phi values (multi circles)
Image Image::initializePhiPart_6(double R1, double R2)
{
	int i, j;
	double x0 = (double)row/4.50;	//cneter of the circle
	double y0 = (double)col/3.50;
	double dx =60.0, dy = 70.0; // dx = 70.0, dy = 90.0; //
	//dy 140 for tri1
	double e =1.;

	double d[3][3];
	double dd[3][3];

	Image temp(row, col);
	allocateMem();


	for (i=0; i<row; i++)
	for (j=0; j<col; j++)
	{
		for (int m=0; m<3; m++)
		for (int n = 0; n<2; n++)
		{
			d[m][n] = sqrt(pow((i - (x0+m*dx)),2)+ pow((j - (y0+n*dy)),2));
			dd[m][n] = sqrt(pow((i - (x0+m*dx)),2)+ pow((j - (y0+n*dy + 15.)),2));
		}															// 50 for tri1

		for (m=0; m<3; m++)
		for (int n = 0; n<2; n++)
		{
			if(d[m][n] > R1-e && d[m][n] <= R1+e)
			temp(i,j)= 155;  // firsrt circle set
			else if (dd[m][n] > R2-e && dd[m][n] <= R2+e)
			temp(i,j)= -255; // second circle set
		}

// initial phi1 and phi2
		for (m=0; m<3; m++)
		for (int n =0; n<2; n++)
		{
			if(d[m][n]<R1-e)
				phi1[i][j]= 1;
			else if (d[m][n] > R1+e)
				phi1[i][j]= -1;
			else phi1[i][j]= 0;

		// initial phi2
			if(dd[m][n]<R2-0)
				phi2[i][j]= 1;
			else if (dd[m][n] > R2+0)
				phi2[i][j]= -1;
			else phi2[i][j]= 0;
		}
	} // end for (i ,j) loop

		return temp;
}

//initial a curve acrroding phi values (multi circles)
Image Image::initializePhiPart_m(double R1, double R2)
{
	int i, j;

/*	double x0 = (double)row/11.0;	//cneter of the circle
	double y0 = (double)col/11.0;
	double dx = 28.0, dy = 36.0;	// 8, 6
*/
	double x0 = (double)row/11.0;	//cneter of the circle
	double y0 = (double)col/11.0;
	double dx = 21.0, dy = 21.0;	// 10, 10

  double e = 0.5;

	double d[10][10];
	double dd[10][10];

	Image temp(row, col);
	allocateMem();


	for (i=0; i<row; i++)
	for (j=0; j<col; j++)
	{
		for (int m=0; m<10; m++)
		for (int n = 0; n<10; n++)
		{
			d[m][n] = sqrt(pow((i - (x0+m*dx)),2)+ pow((j - (y0+n*dy)),2));
			dd[m][n] = sqrt(pow((i - (x0+m*dx)),2)+ pow((j - (y0+n*dy +8.)),2));
		}				//	[8]6], 10

		for (m=0; m<10; m++)
		for (int n = 0; n<10; n++)
		{
			if(d[m][n] > R1-e && d[m][n] <= R1+e)
			temp(i,j)= 155;  // firsrt circle set
			else if (dd[m][n] > R2-e && dd[m][n] <= R2+e)
			temp(i,j)= -255; // second circle set
		}

// initial phi1 and phi2
		int flag=0;
		for (m=0; m<10; m++){
			for (int n =0; n<10; n++){
			if(dd[m][n]<R1){
				phi1[i][j]= 0.1;
				flag=1;
				break;
			}
			}
			if(flag==1)
				break;
		}
		if(flag==0)
			phi1[i][j]= -0.1;
		//	else phi1[i][j]= 0;

		// initial phi2
		flag=0;
		for (m=0; m<10; m++){
			for (int n =0; n<10; n++){
			if(d[m][n]<R2){
				phi2[i][j]= 0.1;
				flag=1;
				break;
			}
			}
			if(flag==1)
				break;
		}
		if(flag==0)
			phi2[i][j]= -0.1;
	} // end for (i ,j) loop

		return temp;
}



//initial a curve acrroding phi values (multi circles)
Image Image::initializePhiPart_n(double R1, double R2)
{
	int i, j;

/*	double x0 = (double)row/11.0;	//cneter of the circle
	double y0 = (double)col/11.0;
	double dx = 28.0, dy = 36.0;	// 8, 6
*/
	double x0 = (double)row/20.0;	//cneter of the circle
	double y0 = (double)col/20.0;
	double dx = 12.0, dy = 12.0;	// 10, 10

  double e = 0.5;

	double d[20][20];
	double dd[20][20];

	Image temp(row, col);
	allocateMem();


	for (i=0; i<row; i++)
	for (j=0; j<col; j++)
	{
		for (int m=0; m<20; m++)
		for (int n = 0; n<20; n++)
		{
			d[m][n] = sqrt(pow((i - (x0+m*dx)),2)+ pow((j - (y0+n*dy)),2));
			dd[m][n] = sqrt(pow((i - (x0+m*dx)),2)+ pow((j - (y0+n*dy +4.)),2));
		}				//	[8]6], 10

		for (m=0; m<20; m++)
		for (int n = 0; n<20; n++)
		{
			if(d[m][n] > R1-e && d[m][n] <= R1+e)
			temp(i,j)= 155;  // firsrt circle set
			else if (dd[m][n] > R2-e && dd[m][n] <= R2+e)
			temp(i,j)= -255; // second circle set
		}

// initial phi1 and phi2
		int flag=0;
		for (m=0; m<20; m++){
			for (int n =0; n<20; n++){
			if(dd[m][n]<R1){
				phi1[i][j]= 0.1;
				flag=1;
				break;
			}
			}
			if(flag==1)
				break;
		}
		if(flag==0)
			phi1[i][j]= -0.1;
		//	else phi1[i][j]= 0;

		// initial phi2
		flag=0;
		for (m=0; m<20; m++){
			for (int n =0; n<20; n++){
			if(d[m][n]<R2){
				phi2[i][j]= 0.1;
				flag=1;
				break;
			}
			}
			if(flag==1)
				break;
		}
		if(flag==0)
			phi2[i][j]= -0.1;
	} // end for (i ,j) loop

		return temp;
}


// calculate the average inside in cirle and outside the cirle
void Image::average(double &n1, double &n2, double &n3, double &n4,
					double &c01, double &c02,double &c03, double &c04 )
{

	double epsilon = 0.01;
	double c1=0., c2=0., c3=0., c4=0. ;
   // n1=0, n2=0, n3=0, n4=0;
//	c01=0., c02=0., c03=0., c04=0.;

	Image Hphi1(row, col), Hphi2(row, col);

	for (int x = 0; x < row; x++)
	for (int y = 0; y < col; y++)
	{
		Hphi1(x,y) = 1.0/2.0*(1+2./PI*atan(phi1[x][y]/epsilon));
		Hphi2(x,y) = 1.0/2.0*(1+2./PI*atan(phi2[x][y]/epsilon));
	}

	for (int i = 0; i < row; i++)
	for (int j = 0; j < col; j++)
	{

		c1 += pixel[i*col+j]*Hphi1(i,j)*Hphi2(i,j);
		n1 += Hphi1(i,j)*Hphi2(i,j);
		c2 += pixel[i*col+j]*(1.-Hphi1(i,j))*(1.-Hphi2(i,j));
		n2 += (1.-Hphi1(i,j))*(1.-Hphi2(i,j));
		c3 += pixel[i*col+j]*Hphi1(i,j)*(1.-Hphi2(i,j));
		n3 += Hphi1(i,j)*(1.-Hphi2(i,j));
		c4 += pixel[i*col+j]*(1.-Hphi1(i,j))*Hphi2(i,j);
		n4 += (1.-Hphi1(i,j))*Hphi2(i,j);

	}

	c01 = c1/n1;
	c02 = c2/n2;
	c03 = c3/n3;
	c04 = c4/n4;
/*
//	int m=0, n=0;
	double c1=0., c2=0., c3=0., c4=0. ;
	n1=0, n2=0, n3=0, n4=0;
	c01=0., c02=0., c03=0., c04=0.;

	for (int i = 0; i < row; i++)
	for (int j = 0; j < col; j++)
	{
		if ((phi1[i][j] > 0)&& (phi2[i][j] > 0)) {		//Inside
			c1 += pixel[i*col + j];
			n1++;
		} else if ((phi1[i][j] < 0)&&(phi2[i][j]< 0)){ //if (phi[i][j] >=0)
			c2 += pixel[i*col + j];
			n2++;
		} else if ((phi1[i][j] > 0) && (phi2[i][j] < 0)) {
			c3 += pixel[i*col + j];
			n3++;
		} else if ((phi1[i][j] < 0)&& (phi2[i][j] > 0)) {
			c4 += pixel[i*col + j];
			n4++;
		}
	}

		c01 = c1/(double)n1;
		c02 = c2/(double)n2;
		c03 = c3/(double)n3;
		c04 = c4/(double)n4;
*/
}



//calculate the difference of energy
double Image::energyDifference(int x, int y, double &m, double &n,
							   double &c1, double &c2)
{
	double F12=0.;
	F12 = pow((pixel[x*col+y]-c2),2)*n/(double)(n+1) -
		  pow((pixel[x*col+y]-c1),2)*m/(double)(m-1);

		return F12;
}

double Image::energyDifference2(int x, int y, double &m, double &n,
							   double &c1, double &c2)
{
	double F21=0.;
	F21 = pow((pixel[x*col+y]-c1),2)*m/(double)(m+1) -
		  pow((pixel[x*col+y]-c2),2)*n/(double)(n-1);

		return F21;
}

// update the image after iterations Using direct energy calculation
Image Image::updateD()
{
	int x, y;

	double n00=0, n10=0, n01=0, n11=0;
	double c00=0., c10=0., c01=0., c11=0.;

	double E11_00=0., E11_10=0., E11_01=0.;
	double E00_11=0., E00_10=0., E00_01=0.;
	double E10_11=0., E10_00=0., E10_01=0.;
	double E01_11=0., E01_10=0., E01_00=0.;

	Image temp1 (row, col);
	Image temp2 (row, col);

	average(n11, n00, n10, n01, c11, c00, c10, c01);

	cout<<"m = " <<n11 <<"	"<<n00  <<"	"<<n10<<"	"<<n01<<endl;
	cout<<"c0 = " <<c11 <<"	"<<c00  <<"	"<<c01<<"	"<<c10<<endl;

	for (x = 0; x < row; x++)
	for (y = 0; y < col; y++)
	{
		if((phi1[x][y] >= 0)&& (phi2[x][y] >= 0))
		{
		 E11_00 = energyDifference(x, y, n11, n00, c11, c00);
		 E11_10 = energyDifference(x, y, n11, n10, c11, c10);
		 E11_01 = energyDifference(x, y, n11, n01, c11, c01);
			double minE =0.;
			minE =min3(E11_00, E11_01, E11_10);
		  if (minE<0) {
			if( E11_00<= minE)
			{
				phi1[x][y] = -1;
				phi2[x][y] = -1;
			}else if( E11_10 <= minE)
			{
				phi1[x][y] = 1;
				phi2[x][y] = -1;
			}else if  (E11_01<= minE)
			{
				phi1[x][y] = -1;
				phi2[x][y] = 1;
			}
		  } else //if  (minE>0)
			{
				phi1[x][y] = 1; //1
				phi2[x][y] = 1;
			}
		} else if ((phi1[x][y] < 0)&&(phi2[x][y] < 0))
		{
		   E00_11 = energyDifference2(x, y, n11, n00, c11, c00);
		   E00_10 = energyDifference(x, y, n00, n10, c00, c10);
		   E00_01 = energyDifference(x, y, n00, n01, c00, c01);

			double minE =0.;
			minE =min3(E00_11, E00_01, E00_10);
		  if (minE <0) {
			if( E00_11<=minE )
			{
				phi1[x][y] = 1;
				phi2[x][y] = 1;
			}else if( E00_10 <= minE)
			{
				phi1[x][y] = 1;
				phi2[x][y] = -1;
			}else if  (E00_01<= minE)
			{
				phi1[x][y] = -1;
				phi2[x][y] = 1;
			}
		  } else //if  (minE>0)
			{
				phi1[x][y] = -1;
				phi2[x][y] = -1;
			}
		} else if ((phi1[x][y] >= 0)&&(phi2[x][y] < 0))
		{
		   E10_11 = energyDifference(x, y, n10, n11, c10, c11);
		   E10_00 = energyDifference(x, y, n10, n00, c10, c00);
		   E10_01 = energyDifference(x, y, n10, n01, c10, c01);
			double minE =0.;
			minE =min3(E10_11, E10_00, E10_01);
		  if (minE <0) {
			if( E10_11<= minE)
			{
				phi1[x][y] = 1;
				phi2[x][y] = 1;
			}else if( E10_00 <= minE)
			{
				phi1[x][y] = -1;
				phi2[x][y] = -1;
			}else if  (E10_01<= minE)
			{
				phi1[x][y] = -1;
				phi2[x][y] = 1;
			}
		  } else //if  (minE>0)
			{
				//if
			    phi1[x][y] =  1;
				phi2[x][y] = -1;
			}
		}else if ((phi1[x][y] < 0)&&(phi2[x][y] >= 0))
		{
		   E01_11 = energyDifference(x, y, n01, n11, c01, c11);
		   E01_00 = energyDifference(x, y, n01, n00, c01, c00);
		   E01_10 = energyDifference(x, y, n01, n10, c01, c10);

			double minE =0.;
			minE =min3(E01_11, E01_10, E01_00);
		  if (minE <0) {
			if( E01_11<= minE)
			{
				phi1[x][y] = 1;
				phi2[x][y] = 1;
			}else if( E01_10 <= minE)
			{
				phi1[x][y] = 1;
				phi2[x][y] = -1;
			}else if  (E01_00<= minE)
			{
				phi1[x][y] = -1;
				phi2[x][y] = -1;
			}
		  } else //if(minE>0)
			{
				phi1[x][y] = -1;
				phi2[x][y] = 1;
			}
		}
	}

	average(n11, n00, n10, n01, c11, c00, c10, c01);

	cout<<"n = " <<n11 <<"	"<<n00  <<"	"<<n10<<"	"<<n01<<endl;
	cout<<"c = " <<c11 <<"	"<<c00  <<"	"<<c10<<"	"<<c01<<endl;

//	Segmentation Only

	for (x = 0; x < row; x++)
	for (y = 0; y < col; y++)
	{
		if((phi1[x][y] >0)&&(phi2[x][y] >0) )
		{
			if((x > 0 && y > 0)&&
			(((phi1[x][y-1] != phi1[x][y]) || (phi2[x][y-1] != phi2[x][y]))
			||((phi1[x-1][y]!= phi1[x][y]) || (phi2[x-1][y] != phi2[x][y]))))
				temp1(x,y)=  c11;//0;//pixel[x*col + y]; // 0;
			else
			temp1(x,y)=	c11;
		} else if((phi1[x][y] <0)&&(phi2[x][y] <0) )
		{
			if((x > 0 && y > 0) &&
			(((phi1[x][y-1] != phi1[x][y]) || (phi2[x][y-1] != phi2[x][y]))
			||((phi1[x-1][y]!= phi1[x][y]) || (phi2[x-1][y] != phi2[x][y]))))
				temp1(x,y)= c00;//c00;//0;//pixel[x*col + y]; //0;
			else
				temp1(x,y)=	c00;//00;
		} else if((phi1[x][y] >0)&&(phi2[x][y] <0) )
		{
			if((x > 0 && y > 0) &&
			(((phi1[x][y-1] != phi1[x][y]) || (phi2[x][y-1] != phi2[x][y]))
			||((phi1[x-1][y]!= phi1[x][y]) || (phi2[x-1][y] != phi2[x][y]))))
				temp1(x,y)= c10;//pixel[x*col + y]; //0;
			else
				temp1(x,y)=	c10;
		} else if((phi1[x][y] <0)&&(phi2[x][y] >0) )
		{
			if((x > 0 && y > 0) &&
			(((phi1[x][y-1] != phi1[x][y]) || (phi2[x][y-1] != phi2[x][y]))
			||((phi1[x-1][y]!= phi1[x][y]) || (phi2[x-1][y] != phi2[x][y]))))
				temp1(x,y)= pixel[x*col + y]; //0;
			else
				temp1(x,y)=	c01;
		}
		//	else if( phi1[x][y] == 0 || phi2[x][y]==0)
		//	temp1(x,y)= 255;
	}


	for (x = 0; x < row-0; x++)
	for (y = 0; y < col-0; y++)
	{
		if ( (x > 0 && y > 0 )&& //x<row-1 && y<col-1) &&
		 ( ((phi1[x][y-1] != phi1[x][y]) || (phi2[x][y-1] != phi2[x][y]))
		 ||((phi1[x-1][y] != phi1[x][y]) || (phi2[x-1][y] != phi2[x][y]))
		 ||((phi1[x-1][y] != phi1[x][y]) || (phi2[x-1][y] != phi2[x][y]))
		 ||((phi1[x][y-1] != phi1[x][y]) || (phi2[x][y-1] != phi2[x][y])) )
		 )
		  temp2(x,y)= 255;
		else
		 temp2(x, y) = temp1(x, y);
	}

//		return temp1; // no contours
		return temp2;
}

// segmentation then denoising
Image Image::DeNoise(Image &temp0)
{
	int x, y;

	double n00=0, n10=0, n01=0, n11=0;
	double c00=0., c10=0., c01=0., c11=0.;

	double E11_00=0., E11_10=0., E11_01=0.;
	double E00_11=0., E00_10=0., E00_01=0.;
	double E10_11=0., E10_00=0., E10_01=0.;
	double E01_11=0., E01_10=0., E01_00=0.;

	Image temp1 (row, col);
	Image temp2 (row, col);

	average(n11, n00, n10, n01, c11, c00, c10, c01);

	cout<<"m = " <<n11 <<"	"<<n00  <<"	"<<n10<<"	"<<n01<<endl;
	cout<<"c0 = " <<c11 <<"	"<<c00  <<"	"<<c01<<"	"<<c10<<endl;

// segementation
	for (x = 0; x < row; x++)
	for (y = 0; y < col; y++)
	{
		if((phi1[x][y] >= 0)&& (phi2[x][y] >= 0))
		{
		 E11_00 = energyDifference(x, y, n11, n00, c11, c00);
		 E11_10 = energyDifference(x, y, n11, n10, c11, c10);
		 E11_01 = energyDifference(x, y, n11, n01, c11, c01);
			double minE =0.;
			minE =min3(E11_00, E11_01, E11_10);
		  if (minE<0) {
			if( E11_00<= minE)
			{
				phi1[x][y] = -1;
				phi2[x][y] = -1;
			}else if( E11_10 <= minE)
			{
				phi1[x][y] = 1;
				phi2[x][y] = -1;
			}else if  (E11_01<= minE)
			{
				phi1[x][y] = -1;
				phi2[x][y] = 1;
			}
		  } else //if  (minE>0)
			{
				phi1[x][y] = 1;
				phi2[x][y] = 1;
			}
		} else if ((phi1[x][y] < 0)&&(phi2[x][y] < 0))
		{
		   E00_11 = energyDifference2(x, y, n11, n00, c11, c00);
		   E00_10 = energyDifference(x, y, n00, n10, c00, c10);
		   E00_01 = energyDifference(x, y, n00, n01, c00, c01);

			double minE =0.;
			minE =min3(E00_11, E00_01, E00_10);
		  if (minE <0) {
			if( E00_11<=minE )
			{
				phi1[x][y] = 1;
				phi2[x][y] = 1;
			}else if( E00_10 <= minE)
			{
				phi1[x][y] = 1;
				phi2[x][y] = -1;
			}else if  (E00_01<= minE)
			{
				phi1[x][y] = -1;
				phi2[x][y] = 1;
			}
		  } else //if  (minE>0)
			{
				phi1[x][y] = -1;
				phi2[x][y] = -1;
			}
		} else if ((phi1[x][y] >= 0)&&(phi2[x][y] < 0))
		{
		   E10_11 = energyDifference(x, y, n10, n11, c10, c11);
		   E10_00 = energyDifference(x, y, n10, n00, c10, c00);
		   E10_01 = energyDifference(x, y, n10, n01, c10, c01);
			double minE =0.;
			minE =min3(E10_11, E10_00, E10_01);
		  if (minE <0) {
			if( E10_11<= minE)
			{
				phi1[x][y] = 1;
				phi2[x][y] = 1;
			}else if( E10_00 <= minE)
			{
				phi1[x][y] = -1;
				phi2[x][y] = -1;
			}else if  (E10_01<= minE)
			{
				phi1[x][y] = -1;
				phi2[x][y] = 1;
			}
		  } else //if  (minE>0)
			{
			    phi1[x][y] =  1;
				phi2[x][y] = -1;
			}
		}else if ((phi1[x][y] < 0)&&(phi2[x][y] >= 0))
		{
		   E01_11 = energyDifference(x, y, n01, n11, c01, c11);
		   E01_00 = energyDifference(x, y, n01, n00, c01, c00);
		   E01_10 = energyDifference(x, y, n01, n10, c01, c10);

			double minE =0.;
			minE =min3(E01_11, E01_10, E01_00);
		  if (minE <0) {
			if( E01_11<= minE)
			{
				phi1[x][y] = 1;
				phi2[x][y] = 1;
			}else if( E01_10 <= minE)
			{
				phi1[x][y] = 1;
				phi2[x][y] = -1;
			}else if  (E01_00<= minE)
			{
				phi1[x][y] = -1;
				phi2[x][y] = -1;
			}
		  } else //if(minE>0)
			{
				phi1[x][y] = -1;
				phi2[x][y] = 1;
			}
		}
	}
// output the segmented 4 regions constants
	average(n11, n00, n10, n01, c11, c00, c10, c01);

 cout<<"n = " <<n11 <<"	"<<n00  <<"	"<<n10<<"	"<<n01<<endl;
 cout<<"c = " <<c11 <<"	"<<c00  <<"	"<<c10<<"	"<<c01<<endl;


//DeNoising


	Image oimg11(row, col), oimg00(row, col), oimg10(row, col),
		  oimg01(row, col); //, oimge(row, col);
	heatADT  Xn11, Xn00, Xn10, Xn01;	//	Xn,

	int nx = row;
	int ny = col;
	double kapp = 1.0;
	double dtt = .50;

	Matrix uRn11(nx, ny);
	Matrix uRnewn11(nx, ny);

	Matrix uRn00(nx, ny);
	Matrix uRnewn00(nx, ny);

	Matrix uRn10(nx, ny);
	Matrix uRnewn10(nx, ny);

	Matrix uRn01(nx, ny);
	Matrix uRnewn01(nx, ny);


	for (x = 0; x < row; x++)
	for (y = 0; y < col; y++)
	{
		if (phi1[x][y] >=0 &&phi2[x][y] >=0)
		{
			if(x>0 && y>0 && x<row-1 && y<col-1)
			{
				if(phi1[x][y-1] != phi1[x][y]||phi2[x][y-1] != phi2[x][y]) // boundary conditions
				{	oimg11(x,y)	= temp0(x,y-1);
				}
				else if(phi1[x-1][y] != phi1[x][y]||phi2[x-1][y] != phi2[x][y])
				{	oimg11(x,y)	= temp0(x-1,y);
				}
				else if(phi1[x][y+1] != phi1[x][y]||phi2[x][y+1] != phi2[x][y])
				{
					oimg11(x,y)	= temp0(x,y+1);
				}
				else if(phi1[x+1][y] != phi1[x][y]||phi2[x+1][y] != phi2[x][y])
				{
					oimg11(x,y)	= temp0(x+1,y);
				} // end of boundary conditions
				else{
				oimg11(x,y)	= temp0(x,y);
				}
			} else
			oimg11(x,y)	= temp0(x,y);

			oimg00(x,y)	= c00;
			oimg10(x,y)	= c10;
			oimg01(x,y)	= c01;
		} else if ((phi1[x][y] < 0)&&(phi2[x][y] < 0))
		{
			oimg00(x,y)	= temp0(x,y);
			oimg11(x,y)	= c11;
			oimg10(x,y)	= c10;
			oimg01(x,y)	= c01;
		}else if ((phi1[x][y] >= 0)&&(phi2[x][y] < 0))
		{
			if(x>0 && y>0 && x<row-1 && y<col-1)
			{
				if(phi1[x][y-1] != phi1[x][y]) // boundary conditions
				{	oimg10(x,y)	= temp0(x,y-1);
				}
				else if(phi1[x-1][y] != phi1[x][y])
				{	oimg10(x,y)	= temp0(x-1,y);
				}
				else if(phi1[x][y+1] != phi1[x][y])
				{
					oimg10(x,y)	= temp0(x,y+1);
				}
				else if(phi1[x+1][y] != phi1[x][y])
				{
					oimg10(x,y)	= temp0(x+1,y);
				} // end of boundary conditions
				else{
				oimg10(x,y)	= temp0(x,y);
				}
			} else
				oimg10(x,y)	= temp0(x,y);

				oimg00(x,y)	= c00;
				oimg11(x,y)	= c11;
				oimg01(x,y)	= c01;
		} else if ((phi1[x][y] < 0)&&(phi2[x][y] >= 0))
		{
			if(x>0 && y>0 && x<row-1 && y<col-1)
			{
				if(phi2[x][y-1] != phi2[x][y]) // boundary conditions
				{	oimg01(x,y)	= temp0(x,y-1);
				}
				else if(phi2[x-1][y] != phi2[x][y])
				{	oimg01(x,y)	= temp0(x-1,y);
				}
				else if(phi2[x][y+1] != phi2[x][y])
				{
					oimg01(x,y)	= temp0(x,y+1);
				}
				else if(phi2[x+1][y] != phi2[x][y])
				{
					oimg01(x,y)	= temp0(x+1,y);
				} // end of boundary conditions
				else{
				oimg01(x,y)	= temp0(x,y);
				}
			} else
				oimg01(x,y)	= temp0(x,y);

			oimg00(x,y)	= c00;
			oimg11(x,y)	= c11;
			oimg10(x,y)	= c10;
			}
	}

	Xn11.init(uRn11, oimg11);
	Xn00.init(uRn00, oimg00);
	Xn10.init(uRn10, oimg10);
	Xn01.init(uRn01, oimg01);

// for phi1 > 0 and phi2 > 0
	for (int k=0; k<3; k++)
	{
		uRnewn11 = Xn11.solver(kapp, dtt, uRn11); // solve Gassian denoise equation using ADI
		uRn11 = uRnewn11;
	}
// for phi1 < 0 and phi2 < 0
	for (k=0; k<2; k++)
	{
		uRnewn00 = Xn00.solver(kapp, dtt, uRn00); // solve Gassian denoise equation using ADI
		uRn00 = uRnewn00;
	}
// for phi1 > 0 and phi2 < 0
	for (k=0; k<2; k++)
	{
		uRnewn10 = Xn10.solver(kapp, dtt, uRn10); // solve Gassian denoise equation using ADI
		uRn10 = uRnewn10;
	}
// for phi1 < 0 and phi2 > 0
	for (k=0; k<3; k++)
	{
		uRnewn01 = Xn01.solver(kapp*10, dtt, uRn01); // solve Gassian denoise equation using ADI
		uRn01 = uRnewn01;
	}

// Output the Denoised image

	for (x = 0; x < row; x++)
	for (y = 0; y < col; y++)
	{
		if((phi1[x][y] >0)&&(phi2[x][y] >0) )
		{
 			temp1(x,y)=	uRn11(x+1,y+1);
		} else if((phi1[x][y] <0)&&(phi2[x][y] <0) )
		{
			temp1(x,y)= uRn00(x+1,y+1);
		} else if((phi1[x][y] >=0)&&(phi2[x][y] <0) )
		{
			temp1(x,y)=	 uRn10(x+1,y+1);
		} else if((phi1[x][y] <0)&&(phi2[x][y] >=0) )
		{
			temp1(x,y)=	uRn01(x+1,y+1);
		}
	}


 // Output and mark the edges

	for (x = 0; x < row-0; x++)
	for (y = 0; y < col-0; y++)
	{
		if ( (x > 0 && y > 0 && x<row-1 && y<col-1) &&
		 ( ((phi1[x][y-1] != phi1[x][y]) || (phi2[x][y-1] != phi2[x][y]))
		 ||((phi1[x-1][y] != phi1[x][y]) || (phi2[x-1][y] != phi2[x][y]))
		 ||((phi1[x+1][y] != phi1[x][y]) || (phi2[x+1][y] != phi2[x][y]))
		 ||((phi1[x][y+1] != phi1[x][y]) || (phi2[x][y+1] != phi2[x][y])))
		 )
		  temp2(x,y)= 255;
		else
		 temp2(x, y) = temp1(x, y);
	}

		return temp1; // no contours oimg01; //
//		return temp2;
}

// ***************************************************************************
//
//	Solving the level set equations
//
// ***************************************************************************


// calculate the average inside in cirle and outside the cirle
void Image::average(double &c11, double &c00,double &c10, double &c01)
{
	double epsilon = 0.01;
	double c1=0., c2=0., c3=0., c4=0. ;
	double n1=0, n2=0, n3=0, n4=0;
	c11=0., c00=0., c10=0., c01=0.;

	Image Hphi1(row, col), Hphi2(row, col);

	for (int x = 0; x < row; x++)
	for (int y = 0; y < col; y++)
	{
		Hphi1(x,y) = 1.0/2.0*(1+2./PI*atan(phi1[x][y]/epsilon));
		Hphi2(x,y) = 1.0/2.0*(1+2./PI*atan(phi2[x][y]/epsilon));
	}

	for (int i = 0; i < row; i++)
	for (int j = 0; j < col; j++)
	{

		c1 += pixel[i*col+j]*Hphi1(i,j)*Hphi2(i,j);
		n1 += Hphi1(i,j)*Hphi2(i,j);
		c2 += pixel[i*col+j]*(1.-Hphi1(i,j))*(1.-Hphi2(i,j));
		n2 += (1.-Hphi1(i,j))*(1.-Hphi2(i,j));
		c3 += pixel[i*col+j]*Hphi1(i,j)*(1.-Hphi2(i,j));
		n3 += Hphi1(i,j)*(1.-Hphi2(i,j));
		c4 += pixel[i*col+j]*(1.-Hphi1(i,j))*Hphi2(i,j);
		n4 += (1.-Hphi1(i,j))*Hphi2(i,j);

	}

	c11 = c1/n1;
	c00 = c2/n2;
	c10 = c3/n3;
	c01 = c4/n4;
}

// calculate the average inside in cirle and outside the cirle
void Image::average(double &c1, double &c0)
{
	double n1=0, n0=0;
	c1=0;
	c0=0;

	for (int i = 0; i < row; i++)
	for (int j = 0; j < col; j++)
	{
		if(phi1[i][j]>=0){
			c1+=pixel[i*col+j];
			n1++;
		}
		else{
			c0+=pixel[i*col+j];
			n0++;
		}
	}

	c1 = c1/n1;
	c0 = c0/n0;
}

// update the image after iterations Using Level set Method
Image Image::updateLevel(double &fact)
{

	int x, y;

// segmentation

	Image temp (row, col);
	Image temp1(row, col);
	Image temp2(row, col);

	LevelADT X1(row,col), X2(row,col);
	int nx = row;
	int ny = col;
	double dtt = 0.005;

	Matrix  Hphi2(row, col), Hphi1(row, col);

	for (x = 0; x < row; x++)
	for (y = 0; y < col; y++)
	{
		temp(x,y) = pixel[x*col+y];
	}


	Matrix uR1(nx, ny);
	Matrix uRnew1(nx, ny);
	Matrix uR2(nx, ny);
	Matrix uRnew2(nx, ny);

	for (x = 0; x < row; x++)
	for (y = 0; y < col; y++)
	{
		uR1(x+1,y+1) = phi1[x][y];
		uR2(x+1,y+1) = phi2[x][y];
	}
	double c11,c00,c10,c01;
	average1(c11, c00, c10, c01);
	double ca1 = 0., cb1 = 0., ca2 = 0., cb2 = 0.,
		   ca3 = 0., cb3 = 0., ca4 = 0., cb4 = 0.;
	int tstep2=1;
//---------------------------------------------
	char WinName[]="Evolution";
/*	char videoName[]="v.avi";

	CvVideoWriter* writer;
	if((writer=cvCreateVideoWriter(videoName,-1,
		10,cvSize(col,row)))==0){
		printf("Output video can't be created.");
		return temp1;
	}
*/
    cvNamedWindow( WinName, CV_WINDOW_AUTOSIZE );

	IplImage* iout = cvCreateImage(cvSize(col,row), 8, 3 );
//------------------------------------------------
	trans(iout);
	cvSaveImage("init.tif",iout);

	int n=0;
	while(tstep2 >=1)
	{
		n++;
		average1(c11, c00, c10, c01);

		for (x = 0; x < row; x++)
		for (y = 0; y < col; y++)
		{
/*			Hphi1(x+1,y+1) = 1.0/2.0*(1+2./PI*atan(phi1[x][y]/epsilon));
			Hphi2(x+1,y+1) = 1.0/2.0*(1+2./PI*atan(phi2[x][y]/epsilon));
*/
			if(phi1[x][y]>=0)
				Hphi1(x+1,y+1)=1;
			else
				Hphi1(x+1,y+1)=0;
			if(phi2[x][y]>=0)
				Hphi2(x+1,y+1)=1;
			else
				Hphi2(x+1,y+1)=0;

		}

		X1.phi_1solver(mu, dtt, uR1, temp, c11, c00, c10, c01, Hphi2);
//		uRnew1= X1.ADTsolver1(mu, dt, uR1, temp, c11, c00, c10, c01, Hphi2);
		X1.phi_2solver(mu, dtt, uR2, temp, c11, c00, c10, c01, Hphi1);
//		uRnew2= X1.ADTsolver2(mu, dt, uR2, temp, c11, c00, c10, c01, Hphi1);

		ca1 = c11, ca2 = c00, ca3 = c10, ca4 = c01;

		for (x = 0; x < row; x++)
		for (y = 0; y < col; y++)
		{
			phi1[x][y] = uR1(x+1,y+1);
			phi2[x][y] = uR2(x+1,y+1);
		}
//-----------------------------------------------
		trans(iout);
		cvShowImage(WinName,iout);
//		cvWriteFrame(writer,iout);
//-----------------------------------------------

		average1(c11, c00, c10, c01);
		cb1 = c11, cb2 = c00, cb3 = c10, cb4 = c01;

		if( (fabs(ca1-cb1)< fact)&& (fabs(ca2-cb2)< fact)&&
			(fabs(ca3-cb3)< fact)&&(fabs(ca4-cb4)< fact) )
//		if(n>50)
			tstep2 = 0;
	} // for step loop

	trans(iout);
	cvSaveImage("result.tif",iout);
		return temp2;
}

void Image::updateLevel1(double &fact)
{
	int x, y;
	double c1=0., c0=0.;
	int imax = (int)fact;
	double epsilon = 1.0;

// segmentation

	Image temp (row, col);

	LevelADT X1(row,col), X2(row,col);
	int nx = row;
	int ny = col;
	double dtt = 0.005;


	for (x = 0; x < row; x++)
	for (y = 0; y < col; y++)
	{
		temp(x,y) = pixel[x*col+y];
	}


	Matrix uR1(nx, ny),uR0(nx,ny);

	for (x = 0; x <row; x++)
	for (y = 0; y <col; y++)
	{
		uR1(x+1,y+1) = phi1[x][y];
	}
	uR0=uR1;
	average(c1, c0);
	int tstep2=1;

//---------------------------------------------
	char WinName[]="Evolution";
/*	char videoName[]="v.avi";

	CvVideoWriter* writer;
	if((writer=cvCreateVideoWriter(videoName,-1,
		10,cvSize(col,row)))==0){
		printf("Output video can't be created.");
		return;
	}
*/
    cvNamedWindow( WinName, CV_WINDOW_AUTOSIZE );

	IplImage* iout = cvCreateImage(cvSize(col,row), 8, 3 );
//------------------------------------------------
	trans1(iout);
	cvSaveImage("init.tif",iout);

	double e0,e1;
	e0=energy1();
	while(tstep2)
	{
		average(c1, c0);

		X1.FTCSsolver(mu,dtt,uR1,temp,c1,c0); //FTC approach


		for (x = 0; x < row; x++)
		for (y = 0; y < col; y++)
		{
			phi1[x][y] = uR1(x+1,y+1);
		}
//-----------------------------------------------
		trans1(iout);
		cvShowImage(WinName,iout);
//		cvWriteFrame(writer,iout);
//-----------------------------------------------
/*		tstep2=0;
		for (x = 1; x <= row; x++){
			for (y = 1; y <= col; y++){
				if(uR1(x,y)-uR0(x,y)>fact){
					tstep2=1;
					break;
				}
			}
			if(tstep2==1)
				break;
		}
*/		uR0=uR1;
		e1=energy1();

		if(fabs(e1-e0)< fact)
		   tstep2 = 0;
		else
			e0=e1;
	} // for step loop

	trans1(iout);
	cvSaveImage("result.tif",iout);
//---------------------------------------
	cvReleaseImage(&iout);
	cvDestroyWindow(WinName);
//	cvReleaseVideoWriter(&writer);
//----------------------------------------

		return;
}


Image Image::creatImage(double R1, double R2)
{
	int i, j;
	double dis1, dis2;
	double e = 30.;
//	double x, y, a, c;

//	initializePhi();
	Image temp(row, col);
//	allocateMem();

	for (i=0; i<row; i++)
	for (j=0; j<col; j++)
	{
//*
		dis1 = sqrt(pow((i - (double)row/2.0),2)+ pow((j - (double)col/2.0),2));
		dis2 = sqrt(pow((i - (double)row/4.0),2)+ pow((j - (double)col/4.0),2));
/*
		if(dis1 > R1 &&dis1 <= R1+e )
			temp(i,j)= 80;
		else if (dis1 > R1+e && dis1 <= R1+ 2*e)
			temp(i,j)= 160;
		else if (dis1 > R1+ 2*e && dis1 <= R1+ 3*e)
			temp(i,j)= 240;
		else
			temp(i,j)= 120;

/*/
		if(sqrt(dis1) > R1-e && sqrt(dis1) <= R1+e )
			temp(i,j)= 180;
		else if (sqrt(dis2) >= 0 && sqrt(dis2) <= R2+e )
			temp(i,j)= 40;
		else if (i>row/2 && i <row/1.2&& j>col/6 && j<col-60)
			temp(i,j)= 250;
		else
			temp(i,j)= 200;
/* /

		if ((j-i >= 0) && (j+i >= 255))
				temp(i,j)= 60;
		else if (i <128)
			temp(i,j)= 140;
		else if (i>=128 && j <128)
			temp(i,j)= 100;
		else
			temp(i,j) = 240;

* /
		if(dis1 <= R1)
		{
		if ((j-i >= 0) && (i>= 128))
			temp(i,j)= 30;
		else if ((j+i>= 255)&&(i <128))
			temp(i,j)= 70;
		else if ((i <128)&&(j>=128))
			temp(i,j)= 170;
		else if ((j-i > 0)&&(j<128))
			temp(i,j)= 170;
		else if ((j-i < 0)&&(i <128))
			temp(i,j)= 110;
		else if ((j+i < 255)&&(i >= 128))
			temp(i,j)= 140;

		else if (i>=128 && j <128)
			temp(i,j)= 240;
		else
			temp(i,j) = 30;
		}else
			temp(i,j) = 200;

/ *
		if(dis1 <= R1)
		{
		if ((j-i >= 0) && (j+i>= 255))
			temp(i,j)= 20;
		else if ((j+i< 255)&&(j>= 128))
			temp(i,j)= 200;
		else if ((j <128)&&(j+i < 255))
			temp(i,j)= 60;
		else if ((j-i < 0)&&(j<128))
			temp(i,j)= 245;
/ *		else if ((j-i < 0)&&(i <128))
			temp(i,j)= 110;
		else if ((j+i < 255)&&(i >= 128))
			temp(i,j)= 140;

		else if (i>=128 && j <128)
			temp(i,j)= 240;
	* /	else
			temp(i,j) = 160;
		}else
			temp(i,j) = 110;
/*

		if (i>10&&j>10&&i<240&&j<240)
		{
			if ((j-i >= 0) && (j+i< 250) && i< 158)
				temp(i,j)= 0;
			else if (i>160 &&i<200 && j >40 && j<210)
				 temp(i,j)= 255;
			else
			 temp(i,j)= 140;
		}
		else
			temp(i,j) = 140;
*/
	}

		return temp;
}

void Image::trans(IplImage*iout){
	for(int y=0;y<row;y++){
		for (int x=0;x<col;x++){
			if ((x>1 && y>1 &&  x<row-2 && y<col-2) &&
				(phi1[y][x-1]<0 ||phi1[y-1][x]<0 ||
				phi1[y][x+1]<0 ||phi1[y+1][x]<0)&& phi1[y][x]>0){
				((uchar*)(iout->imageData
					+iout->widthStep*y))[x*3]=0;
				((uchar*)(iout->imageData
					+iout->widthStep*y))[x*3+1]=255;
				((uchar*)(iout->imageData
					+iout->widthStep*y))[x*3+2]=0;
			}
			else if ((x>1 && y>1 &&  x<row-2 && y<col-2) &&
				(phi2[y][x-1]<0 ||phi2[y-1][x]<0 ||
				phi2[y][x+1]<0 ||phi2[y+1][x]<0)&& phi2[y][x]>0){
				((uchar*)(iout->imageData
					+iout->widthStep*y))[x*3]=0;
				((uchar*)(iout->imageData
					+iout->widthStep*y))[x*3+1]=0;
				((uchar*)(iout->imageData
					+iout->widthStep*y))[x*3+2]=255;
			}
			else{
				((uchar*)(iout->imageData+iout->widthStep*y))[x*3]=
					((uchar*)(iout->imageData+iout->widthStep*y))[x*3+1]=
					((uchar*)(iout->imageData+iout->widthStep*y))[x*3+2]=pixel[y*col+x];
			}
		}
	}
}

void Image::trans1(IplImage*iout){
	for(int y=0;y<row;y++){
		for (int x=0;x<col;x++){
			if ((x>1 && y>1 &&  x<row-2 && y<col-2) &&
				(phi1[y][x-1]<0 ||phi1[y-1][x]<0 ||
				phi1[y][x+1]<0 ||phi1[y+1][x]<0)&& phi1[y][x]>0){
				((uchar*)(iout->imageData
					+iout->widthStep*y))[x*3]=0;
				((uchar*)(iout->imageData
					+iout->widthStep*y))[x*3+1]=255;
				((uchar*)(iout->imageData
					+iout->widthStep*y))[x*3+2]=0;
			}
			else{
				((uchar*)(iout->imageData+iout->widthStep*y))[x*3]=
					((uchar*)(iout->imageData+iout->widthStep*y))[x*3+1]=
					((uchar*)(iout->imageData+iout->widthStep*y))[x*3+2]=pixel[y*col+x];
			}
		}
	}
}

void Image::average1(double &c11, double &c00,double &c10, double &c01)
{
	double epsilon = 0.01;
	double c1=0., c2=0., c3=0., c4=0. ;
	double n1=0, n2=0, n3=0, n4=0;
	c11=0., c00=0., c10=0., c01=0.;

	Image Hphi1(row, col), Hphi2(row, col);

	for (int x = 0; x < row; x++)
	for (int y = 0; y < col; y++)
	{
//		Hphi1(x,y) = 1.0/2.0*(1+2./PI*atan(phi1[x][y]/epsilon));
//		Hphi2(x,y) = 1.0/2.0*(1+2./PI*atan(phi2[x][y]/epsilon));

		if(phi1[x][y]>=0)
			Hphi1(x,y)=1;
		else
			Hphi1(x,y)=0;
		if(phi2[x][y]>=0)
			Hphi2(x,y)=1;
		else
			Hphi2(x,y)=0;

	}

	for (int i = 0; i < row; i++)
	for (int j = 0; j < col; j++)
	{

		c1 += pixel[i*col+j]*Hphi1(i,j)*Hphi2(i,j);
		n1 += Hphi1(i,j)*Hphi2(i,j);
		c2 += pixel[i*col+j]*(1.-Hphi1(i,j))*(1.-Hphi2(i,j));
		n2 += (1.-Hphi1(i,j))*(1.-Hphi2(i,j));
		c3 += pixel[i*col+j]*Hphi1(i,j)*(1.-Hphi2(i,j));
		n3 += Hphi1(i,j)*(1.-Hphi2(i,j));
		c4 += pixel[i*col+j]*(1.-Hphi1(i,j))*Hphi2(i,j);
		n4 += (1.-Hphi1(i,j))*Hphi2(i,j);

	}
	if(n1!=0)
		c11 = c1/n1;
	else
		c11=0;
	if(n2!=0)
		c00 = c2/n2;
	else
		c00=0;
	if(n3!=0)
		c10 = c3/n3;
	else
		c10=0;
	if(n4!=0)
		c01 = c4/n4;
	else
		c01=0;
}
double Image::energy1(){
    int x,y;

	double U1,U2;

	double sum1=0,sum2=0,l=0;
	int n1=0,n2=0;
	for(y=0;y<row;y++){
		for(x=0;x<col;x++){
			if(phi1[y][x]==1){
				sum1+=pixel[y*col+x];
				n1++;
			}
			else{
				sum2+=pixel[y*col+x];
				n2++;
			}
		}
	}
	U1=(float)sum1/(float)n1;
	U2=(float)sum2/(float)n2;

	sum1=0;
	sum2=0;
	for(y=0;y<row;y++){
		for(x=0;x<col;x++){
			if(phi1[y][x]==1)
				sum1+=pow((pixel[y*col+x]-U1),2);
			else
				sum2+=pow((pixel[y*col+x]-U2),2);
		}
	}
	for(y=1;y<row-1;y++){
		for(x=1;x<col-1;x++){
				if ((phi1[y][x-1]<0 ||phi1[y-1][x]<0 ||
				phi1[y][x+1]<0 ||phi1[y+1][x]<0)&& phi1[y][x]>0)
				l++;
		}
	}
	double v=sum1+sum2+(double)mu*l;
	return v;
}
void Image::smooth(){
	float kernel[11][11],kernel1[11][11];
	int x,y;
	IplImage* out=cvCreateImage(cvSize(col,row),8,1);
	for(y=0;y<11;y++){
		for(x=0;x<11;x++){
			kernel[x][y]=0.5/3.1415926/a1/a1*exp(-0.5*
					(((float)y-5)*((float)y-5)+((float)x-5)*((float)x-5))/a1/a1);
		}
	}
	for(y=0;y<row;y++){
		for(x=0;x<col;x++){
			float p=0,k=0;
			for(int y1=-5;y1<=5;y1++){
				for(int x1=-5;x1<=5;x1++){
					if((y1+y)<1||(y1+y)>=row-1||(x+x1)<1||(x+x1)>=col-1
						||phi1[y+y1][x+x1+1]*phi1[y][x]<0
						||phi2[y+y1][x+x1+1]*phi2[y][x]<0)
						p+=kernel[5+y1][5+x1]*pixel[y*col+x];
					else
						p+=kernel[5+y1][5+x1]*pixel[(y+y1)*col+x+x1];
					k+=kernel[5+y1][5+x1];
				}
			}
			p/=k;
			if(p<=255)
				((unsigned char*)(out->imageData+out->widthStep*y))[x]=(unsigned char)p;
			else
				((unsigned char*)(out->imageData+out->widthStep*y))[x]=255;
		}
	}
	cvSaveImage("smooth.tif",out);
	cvReleaseImage(&out);
}
