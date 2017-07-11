#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// %     ************
// %
// %     Description
// %     ----------
// %     Compute the Time Warp Edit distance between time series ts1 and ts2.
// % 
// %     Parameters
// %     ----------
// %     param ts1   :  n1 x dim, time series 1 matrix with the length of n1
// %     param ts2   :  n2 x dim, time series 2 matrix with the length of n2
// %
// %     Options
// %     ----------
// %     opt lambda    : double, penalty, punishment for distances at deletions
// %     opt nu        : double, stiffness, determines the elasticity Nu > = 0 required for distance measurement.
// %     opt winSize   : integer, temporal constraint on the warping window
// %                              size. default value = -1
// %
// %     Returns
// %     -------
// %     twedDist      : double, The Time Warp Edit distance between time series ts1 and ts2 
// %
// %     Other m-files required    : TS_euclideanDistance
// %     Subfunctions              : none
// %     MAT-files required        : none
// %     
// %     References
// %     ----------
// %     http://people.irisa.fr/Pierre-Francois.Marteau/
// %     @article{marteau2009time,
// %       title={Time warp edit distance with stiffness adjustment for time series matching},
// %       author={Marteau, Pierre-Fran{\c{c}}ois},
// %       journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
// %       volume={31},
// %       number={2},
// %       pages={306--318},
// %       year={2009},
// %       publisher={IEEE}
// %     }
// %
// %     Author
// %     ----------
// %     Amir Salarpour, Ph.D. Candidate, Artificial Intelligence
// %     Bu-Ali Sina University, Hamedan, Iran, Dept. of Computer Engineering
// %     email address : amir.salarpour@gmail.com  
// %     Website       : http://www.salarpour.com
// %     December 2016 : Last revision: 26-Jan-2017
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

void twedCalculate(double *ts1, int n1, double *ts2, int n2, double lambda, double nu, int winSize, int dim, double *twedDist)
{
	double D[n1 + 1][n2 + 1];
	double Z[dim];
	int i, j;
	int jS, jF;
    int lenDiff;
	double tmp, tmp1, tmpMin;
	
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
	
	for ( i = 0; i < dim; i++)
	{
		Z[i] = 0;
	}
	
	for ( i = 0; i <= n1; i++)
	{
		for ( j = 0; j <= n2; j++)
		{
			D[i][j] = 0;
			
			if (i == 0)
			{
				D[i][j] = INFINITY;
			}
			if (j == 0)
			{
				D[i][j] = INFINITY;
			}
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
			if (i > 0 && j > 0)
			{
				tmp = euclideanCalculate(ts1, n1, ts2, n2, dim, i-1, j-1);
				tmp1 = euclideanCalculate(ts1, n1, ts2, n2, dim, i - 2, j - 2);
				D[i][j] = D[i-1][j-1] + tmp + tmp1 + (nu * 2 * abs(i-j));
			}
			else
			{
				tmp = euclideanCalculate(ts1, n1, ts2, n2, dim, i-1, j-1);
				D[i][j] = D[i-1][j-1] +  tmp + (nu * abs(i-j));
			}
			
			if (i > 0)
			{
				tmp = euclideanCalculate(ts1, n1, ts1, n1, dim, i-1, i-2);
				tmpMin = D[i-1][j] + nu + lambda + tmp;
			}
			else
			{
				tmp = euclideanCalculate(ts1, n1, Z, 1, dim, i-1, 0);
				tmpMin = D[i-1][j] + nu + lambda + tmp;
			}
			if ( tmpMin < D[i][j] )
			{
				D[i][j] = tmpMin;
			}
			
			if (j > 0)
			{
				tmp = euclideanCalculate(ts2, n2, ts2, n2, dim, j-1, j-2);
				tmpMin = D[i][j-1] + nu + lambda + tmp;
			}
			else
			{
				tmp = euclideanCalculate(ts2, n2, Z, 1, dim, j-1, 0);
				tmpMin = D[i][j-1] + nu + lambda + tmp;
			}
			if ( tmpMin < D[i][j] )
			{
				D[i][j] = tmpMin;
			}

		}
	}
	

    twedDist[0] = D[n1][n2];
	
	return;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Macros for the ouput and input arguments */
	#define out1 		plhs[0]
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
	else if(nlhs != 1)
	{
		mexErrMsgTxt("Too many output arguments.");
	}
		
	/* Define local variables */
	double *ts1, *ts2;
	double lambda, nu;
	int winSize;
	int n1, dim1, n2, dim2;
	double *twedDist;
	
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
	
	/* check to make sure winSize is a scalar */
	if(nrhs == 2)
    {
		lambda = 1;
		nu = 0;
        winSize = -1;
    }
	else if(nrhs == 3)
    {
		lambda = mxGetScalar(in3);
		nu = 0;
        winSize = -1;

    }
	else if(nrhs == 4)
    {
		lambda = mxGetScalar(in3);
		nu = mxGetScalar(in4);
        winSize = -1;
	}
	else if(nrhs == 5)
    {
		lambda = mxGetScalar(in3);
		nu = mxGetScalar(in4);
        winSize = (int) mxGetScalar(in5);
	}

	

	
	out1 = mxCreateDoubleMatrix( 1, 1, mxREAL);
	
	twedDist = mxGetPr(out1);
	
	twedCalculate(ts1, n1, ts2, n2, lambda, nu, winSize, dim1, twedDist);
    
    
    return;
    
}



