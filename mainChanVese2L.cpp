/*******************************************************************
 * file:    mainChanVese2L.cpp
 * version: 1.0
 * project:  Project of IEEE paper 
 * author:  Song Gao
 * Modification date:    Feb. 10, 2004
 * description: The main program
 *				Level set Methods and Image segmentations
 *				
 *				Chan and Vese Method 
 *				TWO (2) Level set Functions are used.
 *
 *		Direct energy computation method and
		The level set method are implemented in this project.
 *
 *********************************************************************/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iomanip> 
#include <time.h> 
#include "Image.h"

#include <windows.h>

float mu=10;

double getCPUTime() // return process CPU time in seconds
{
	long cre, ext, ker, usr;
	static HANDLE me=GetCurrentProcess();
	GetProcessTimes(me,(LPFILETIME)&cre,(LPFILETIME)&
					ext,(LPFILETIME)&ker,(LPFILETIME)&usr);
	return (ker+usr)*1.0e-7;
}

int main()
{
    char *outfile1;
	char *outfile2;

	double fact = 0.001;

	
	char infile1[50];
	printf("Input image: ");
	scanf("%s",infile1);
	printf("mu=");
	scanf("%f",&mu);

	Image P1, P2;

	P1.readPGMRAW(infile1); // input image  

	P2 = P1 + P1.initializePhiPart_2(80, 80);//
	
	P1.updateLevel1(fact); // 1 level set function 

//	P1.updateLevel(fact); // 2 level set functions
	return 0;
}

