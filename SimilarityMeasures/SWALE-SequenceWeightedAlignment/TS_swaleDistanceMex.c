#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// %     ************
// %
// %     Description
// %     ----------
// %     Compute the Seuence Weighted Alignment distance between time series ts1 and ts2.
// % 
// %     Parameters
// %     ----------
// %     param ts1   :  n1 x dim, time series 1 matrix with the length of n1
// %     param ts2   :  n2 x dim, time series 2 matrix with the length of n2
// %
// %     Options
// %     ----------
// %     opt eps       : double, a threshold for considering distance
// %     opt gapC      : double, gap cost
// %     opt reward    : double, match reward
// %     opt winSize   : integer, temporal constraint on the warping window
// %                              size. default value = -1
// %
// %     Returns
// %     -------
// %     swaleDist       : double, The Seuence Weighted Alignment distance between time series ts1 and ts2
// %
// %     Other m-files required    : TS_euclideanDistance
// %     Subfunctions              : none
// %     MAT-files required        : none
// %     
// %     References
// %     ----------
// %     @inproceedings{morse2007efficient,
// %       title={An efficient and accurate method for evaluating time series similarity},
// %       author={Morse, Michael D and Patel, Jignesh M},
// %       booktitle={Proceedings of the 2007 ACM SIGMOD international conference on Management of data},
// %       pages={569--580},
// %       year={2007},
// %       organization={ACM}
// %     }   
// %
// %     Author
// %     ----------
// %     Amir Salarpour, Ph.D. Candidate, Artificial Intelligence
// %     Bu-Ali Sina University, Hamedan, Iran, Dept. of Computer Engineering
// %     email address : amir.salarpour@gmail.com  
// %     Website       : http://www.salarpour.com
// %     December 2016 : Last revision: 31-Jan-2017
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

void swaleCalculate(double *ts1, int n1, double *ts2, int n2, double eps, double gapC, double reward, int winSize, int dim, double *swaleDist)
{
	double D[n1 + 1][n2 + 1];
	int i, j;
	int jS, jF;
    double euclDist;
    int lenDiff;
	double tmp;
	
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
	for ( i = 0; i <= n1; i++)
	{
		for ( j = 0; j <= n2; j++)
		{
			D[i][j] = 0;
			
			if (i == 0)
			{
				D[i][j] = gapC * j;
			}
			if (j == 0)
			{
				D[i][j] = gapC * i;
			}
		}
	}

	
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
			tmp = euclideanCalculate(ts1, n1, ts2, n2, dim, i-1, j-1);
			
			if (tmp < eps)
			{
				D[i][j] = D[i-1][j-1] + reward;
			}
			else
			{
				D[i][j] = D[i][j-1] + gapC;
				
				if (D[i][j] < D[i-1][j] + gapC)
				{
					D[i][j] = D[i-1][j] + gapC;
				}
			}

			
		}
	}
	

    swaleDist[0] = D[n1][n2];
	
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
	#define in6 		prhs[5]
	
	/* Check the number of arguments */
	if(nrhs < 2 && nrhs > 6)
	{
		mexErrMsgTxt("Wrong number of input arguments.");
	}
	else if(nlhs != 1 && nlhs != 2)
	{
		mexErrMsgTxt("Too many output arguments.");
	}
		
	/* Define local variables */
	double *ts1, *ts2;
	double eps, gapC, reward;
	int winSize;
	int n1, dim1, n2, dim2;
	double *swaleDist;
	
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
		eps = INFINITY;
		gapC = 0;
		reward = 0;
        winSize = -1;
    }
	else if(nrhs == 3)
    {
		eps = mxGetScalar(in3);
		gapC = 0;
		reward = 0;
        winSize = -1;

    }
	else if(nrhs == 4)
    {
		eps = mxGetScalar(in3);
		gapC = mxGetScalar(in4);
		reward = 0;
        winSize = -1;
	}
	else if(nrhs == 5)
    {
		eps = mxGetScalar(in3);
		gapC = mxGetScalar(in4);
		reward = mxGetScalar(in5);
        winSize = -1;
	}
		else if(nrhs == 6)
    {
		eps = mxGetScalar(in3);
		gapC = mxGetScalar(in4);
		reward = mxGetScalar(in5);
        winSize = (int) mxGetScalar(in6);
	}

	

	
	out1 = mxCreateDoubleMatrix( 1, 1, mxREAL);
	out2 = mxCreateDoubleMatrix( 1, 1, mxREAL);
	
	swaleDist = mxGetPr(out1);
	
	swaleCalculate(ts1, n1, ts2, n2, eps, gapC, reward, winSize, dim1, swaleDist);
    
    
    return;
    
}



