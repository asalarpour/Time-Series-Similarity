#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// %     ************
// %
// %     Description
// %     ----------
// %     Compute the Dynamic-Time Warping distance between time series ts1 and ts2.
// % 
// %     Parameters
// %     ----------
// %     param ts1   :   n1 x dim, time series 1 matrix with the length of n1
// %     param ts2   :   n2 x dim, time series 2 matrix with the length of n
// %
// %     Options
// %     ----------
// %     opt winSize   :   integer,   temporal constraint on the warping window
// %                              size. default value = -1 
// %
// %     Returns
// %     -------
// %     dtwDist       :   double, The Dynamic-Time Warping distance between time series ts1 and ts2.
// %     dtwLen        :   integer, Show the length of warping path between time series
// %
// %     Other files required    : none
// %     Subfunctions              : none
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

void dtwCalculate(double *ts1, int n1, double *ts2, int n2, int winSize, int dim, double *dtwDist, double *dtwLen)
{
	double D[n1 + 1][n2 + 1];
	int L[n1 + 1][n2 + 1];
	int i, j;
	int jS, jF;
    double euclDist;
    int lenDiff;
	
	if (n1 - n2 > 0)
	{
		int lenDiff = n1 - n2;
	}
	else
	{
		int lenDiff = n2 - n1;
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
			D[i][j] += euclDist;
			L[i][j] += 1;
			
		}
	}
	
    
    dtwDist[0] = D[n1][n2];
	dtwLen[0] = L[n1][n2];
	
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
	
	/* Check the number of arguments */
	if(nrhs != 2 && nrhs != 3)
	{
		mexErrMsgTxt("Wrong number of input arguments.");
	}
	else if(nlhs != 1 && nlhs != 2)
	{
		mexErrMsgTxt("Too many output arguments.");
	}
		
	/* Define local variables */
	double *ts1, *ts2;
	int winSize;
	int n1, dim1, n2, dim2;
	double *dtwDist, *dtwLen;
	
	/* check to make sure winSize is a scalar */
	if(nrhs == 2)
    {
        winSize = -1;
    }
	else if(nrhs == 3)
    {
        if( !mxIsDouble(in3) || mxIsComplex(in3) ||
                mxGetN(in3) * mxGetM(in3)!=1 )
        {
            mexErrMsgTxt("windows size should be an scalar");
        }

        winSize = (int) mxGetScalar(in3);
    }
	
	ts1 = mxGetPr(in1);
	ts2 = mxGetPr(in2);
	
	n1 = mxGetM(in1);
	dim1 = mxGetN(in1);
	
	n2 = mxGetM(in2);
	dim2 = mxGetN(in2);
	
	if (dim1 != dim2)
	{
		mexErrMsgTxt("Two trajectories dimension must be the same");
	}
	
	out1 = mxCreateDoubleMatrix( 1, 1, mxREAL);
	out2 = mxCreateDoubleMatrix( 1, 1, mxREAL);
	
	dtwDist = mxGetPr(out1);
	dtwLen = mxGetPr(out2);
	dtwCalculate(ts1, n1, ts2, n2, winSize, dim1, dtwDist, dtwLen);
    
    
    return;
    
}



