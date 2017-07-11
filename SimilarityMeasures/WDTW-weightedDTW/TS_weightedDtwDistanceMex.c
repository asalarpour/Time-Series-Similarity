#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// %     ************
// %
// %     Description
// %     ----------
// %     Compute the Weighted Dynamic-Time Warping distance between time series ts1 and ts2.
// % 
// %     Parameters
// %     ----------
// %     param ts1   :  m x dim, time series 1 matrix with the length of m
// %     param ts2   :  n x dim, time series 2 matrix with the length of n
// %
// %     Options
// %     ----------
// %     opt winSize   :   integer,   temporal constraint on the warping window
// %                              size. default value = -1
// %     opt g         :   double,  empirical constant that controls the curvature (slope) of the function
// %                              g = 0 it is constant, g = 0.05 it is linear, g = 0.25 it is sigmoid, g = 3 two distinct weight
// %                              the range of optimal g is between 0.01 to 0.6
// %     opt weightMax : double,  maximum weight to normalze
// %
// %     Returns
// %     -------
// %     wdtwDist       : double, The Weighted Dynamic-Time Warping distance between time series ts1 and ts2
// %     wdtwLen        :   integer, Show the length of warping path between time series
// %
// %     Other m-files required    : TS_euclideanDistance
// %     Subfunctions              : TS_weightFun
// %     MAT-files required        : none
// %     
// %     References
// %     ----------
// %
// %     Author
// %     ----------
// %     Amir Salarpour, Ph.D. Candidate, Artificial Intelligence
// %     Bu-Ali Sina University, Hamedan, Iran, Dept. of Computer Engineering
// %     email address : amir.salarpour@gmail.com  
// %     Website       : http://www.salarpour.com
// %     December 2016 : Last revision: 27-Jan-2017
// %     
// %     ************

double euclideanCalculate(double *ts1, int n1, double *ts2, int n2, int dim, int idx1, int idx2)
{
	double euclDist = 0;
	int i;
	double tmp;
	
	
	for (i = 0; i < dim; i++)
	{
		tmp = ts1[idx1 + n1 * i] - ts2[idx2 + n2 * i];
		euclDist += tmp * tmp;
	}
	euclDist = sqrt(euclDist);
	return euclDist;
}

double weightCalculate(int i, int j, int n1, int n2, double g, double weightMax)
{
	double output;
	double diff, len;
	
	diff = (double) abs(i - j);
	len = (double) (n1 + n2);
	len = len / 2;
	output = weightMax / ( 1 + exp(-g * (diff - (len / 2) ) ) );
/* 	if (output > 0.4 && output <0.6)
	{
		mexPrintf("%f --- %f ---- \n", diff, len);
		mexPrintf("%d --- %d ---- %f ---- \n", i, j, output);
	} */
	
	return output;
}

void wdtwCalculate(double *ts1, int n1, double *ts2, int n2, double g, double weightMax, int winSize, int dim, double *wdtwDist, double *wdtwLen)

{
	double D[n1 + 1][n2 + 1];
	int L[n1 + 1][n2 + 1];
	int i, j;
	int jS, jF;
    double euclDist, weightVal;
	int lenDiff;
	
	if (n1 - n2 > 0)
	{
		lenDiff = n1 - n2;
	}
	else
	{
		lenDiff = n2 - n1;
	}
	
	if (winSize != -1 && winSize < lenDiff)
	{
		winSize = lenDiff;
	}
	
	for (i = 0; i <= n1; i++)
	{
		for (j = 0; j <= n2; j++)
		{
			D[i][j] = INFINITY;
			L[i][j] = 0;
		}
	}
	
	D[0][0] = 0;
	
	for (i = 1; i <= n1; i++)
	{
		if (winSize == -1)
		{
			jS = 1;
			jF = n2;
		}
		else
		{
			jS = i - winSize > 1 ? i - winSize : 1;
			jF = i + winSize < n2 ? i + winSize : n2;
		}
		for (j = jS; j <= jF; j++)
		{
			weightVal = weightCalculate(i, j, n1, n2, g, weightMax);
/* 			if (weightVal > 0.1 && weightVal <0.9)
			{
				mexPrintf("%d --- %d ---- %f ---- \n", i, j, weightVal);
			} */
			euclDist = euclideanCalculate(ts1, n1, ts2, n2, dim, i-1, j-1);
			D[i][j] = D[i-1][j-1];
			L[i][j] = L[i-1][j-1];
			if (D[i][j] > D[i-1][j])
			{
				D[i][j] = D[i-1][j];
				L[i][j] = L[i-1][j];
			}
			if (D[i][j] > D[i][j-1])
			{
				D[i][j] = D[i][j-1];
				L[i][j] = L[i][j-1];
			}
			D[i][j] += weightVal * euclDist;
			L[i][j] += 1;
			
		}
	}
	
    
    wdtwDist[0] = D[n1][n2];
	wdtwLen[0] = L[n1][n2];
	
	return;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Macros for the ouput and input arguments */
	#define out1 		plhs[0]
	#define out2 		plhs[1]
	#define in1 		prhs[0]
	#define in2 		prhs[1]
	#define in3 		prhs[2]
	#define in4 		prhs[3]
	#define in5 		prhs[4]
	
	/* Check the number of arguments */
	if(nrhs < 2 && nrhs > 5)
	{
		mexErrMsgTxt("Wrong number of input arguments.");
	}
	else if(nlhs != 1 && nlhs != 2)
	{
		mexErrMsgTxt("Too many output arguments.");
	}
		
	/* Define local variables */
	double *ts1, *ts2;
	double g, weightMax;
	int winSize;
	int n1, dim1, n2, dim2;
	double *wdtwDist, *wdtwLen;
	
	/* check to make sure winSize is a scalar */
	if(nrhs == 2)
    {
		g = 0.25;
		weightMax = 1;
        winSize = -1;
    }
	else if (nrhs == 3)
	{
		g = mxGetScalar(in3);
		weightMax = 1;
        winSize = -1;
	}
	else if (nrhs == 4)
	{
		g = mxGetScalar(in3);
		weightMax = mxGetScalar(in4);
        winSize = -1;
	}	
	else if(nrhs == 5)
    {
        g = mxGetScalar(in3);
		weightMax = mxGetScalar(in4);

        winSize = (int) mxGetScalar(in5);
    }
	
	ts1 = mxGetPr(in1);
	ts2 = mxGetPr(in2);
	
	n1 = mxGetM(in1);
	dim1 = mxGetN(in1);
	
	n2 = mxGetM(in2);
	dim2 = mxGetN(in2);
	
	if (dim1 != dim2)
	{
		mexErrMsgTxt("Two time series dimension must be the same");
	}
	
	out1 = mxCreateDoubleMatrix( 1, 1, mxREAL);
	out2 = mxCreateDoubleMatrix( 1, 1, mxREAL);
	
	wdtwDist = mxGetPr(out1);
	wdtwLen = mxGetPr(out2);
	wdtwCalculate(ts1, n1, ts2, n2, g, weightMax, winSize, dim1, wdtwDist, wdtwLen);
    
    
    return;
    
}



