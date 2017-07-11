#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// %     ************
// %
// %     Description
// %     ----------
// %     Compute the Edit Distance on Real sequence between time series ts1 and ts2.
// % 
// %     Parameters
// %     ----------
// %     param ts1   :  m x dim, time series 1 matrix with the length of m
// %     param ts2   :  n x dim, time series 2 matrix with the length of n
// %
// %     Options
// %     ----------
// %     opt eps       : double, a threshold for considering distance
// %                             default value = Inf
// %     opt winSize   : integer, temporal constraint on the warping window
// %                              size. default value = -1
// %
// %     Returns
// %     -------
// %     edrDist       : double, The Edit Distance on Real sequence between time series ts1 and ts2.
// %
// %     Other m-files required    : TS_euclideanDistance
// %     Subfunctions              : none
// %     MAT-files required        : none
// %     
// %     References
// %     ----------
// %     @inproceedings{chen2005robust,
// %       title={Robust and fast similarity search for moving object trajectories},
// %       author={Chen, Lei and {\"O}zsu, M Tamer and Oria, Vincent},
// %       booktitle={Proceedings of the 2005 ACM SIGMOD international conference on Management of data},
// %       pages={491--502},
// %       year={2005},
// %       organization={ACM}
// %     }
// %
// %     Author
// %     ----------
// %     Amir Salarpour, Ph.D. Candidate, Artificial Intelligence
// %     Bu-Ali Sina University, Hamedan, Iran, Dept. of Computer Engineering
// %     email address : amir.salarpour@gmail.com  
// %     Website       : http://www.salarpour.com
// %     December 2016 : Last revision: 28-Jan-2017
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

void edrCalculate(double *ts1, int n1, double *ts2, int n2, double eps, int winSize, int dim, double *edrDist)
{
	double D[n1 + 1][n2 + 1];
	int i, j;
	int jS, jF;
    double euclDist;
    int lenDiff;
	double tmp, subcost;
	
	if (n1 - n2 > 0)
	{
		int lenDiff = n1 - n2;
		tmp = (double) n1;
	}
	else
	{
		int lenDiff = n2 - n1;
		tmp = (double) n2;
	}
	
	if (winSize != -1 && winSize < lenDiff)
	{
		winSize = lenDiff;
	}
	
	for (i = 0; i <= n1; i++)
	{
		for (j = 0; j <= n2; j++)
		{
			D[i][j] = 0;
			
			if (i == 0)
			{
				D[i][j] = (double) j;
			}
			if (j == 0)
			{
				D[i][j] = (double) i;
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
			euclDist = euclideanCalculate(ts1, n1, ts2, n2, dim, i-1, j-1);
			if (euclDist < eps)
			{
				subcost = 0;
			}
			else
			{
				subcost = 1;
			}
			D[i][j] = D[i-1][j-1] + subcost;
			
			if (D[i][j] > D[i-1][j] + 1)
			{
				D[i][j] = D[i-1][j] + 1;
			}
			if (D[i][j] > D[i][j-1] + 1)
			{
				D[i][j] = D[i][j-1] + 1;
			}

			
		}
	}
	

    edrDist[0] = D[n1][n2];
	
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
	
	/* Check the number of arguments */
	if(nrhs != 2 && nrhs != 3 && nrhs != 4)
	{
		mexErrMsgTxt("Wrong number of input arguments.");
	}
	else if(nlhs != 1)
	{
		mexErrMsgTxt("Too many output arguments.");
	}
		
	/* Define local variables */
	double *ts1, *ts2;
	double eps;
	int winSize;
	int n1, dim1, n2, dim2;
	double *edrDist;
	
	/* check to make sure winSize is a scalar */
	if(nrhs == 2)
    {
		eps = INFINITY;
        winSize = -1;
    }
	else if(nrhs == 3)
    {
        if( !mxIsDouble(in3) || mxIsComplex(in3) ||
                mxGetN(in3) * mxGetM(in3)!=1 )
        {
            mexErrMsgTxt("threshold value should be an scalar");
        }

        eps = mxGetScalar(in3);
		winSize = -1;
    }
	else if(nrhs == 4)
    {
        if( !mxIsDouble(in3) || mxIsComplex(in3) ||
                mxGetN(in3) * mxGetM(in3)!=1 )
        {
            mexErrMsgTxt("threshold value should be an scalar");
        }

        eps = mxGetScalar(in3);
		
		if( !mxIsDouble(in4) || mxIsComplex(in4) ||
                mxGetN(in4) * mxGetM(in4)!=1 )
        {
            mexErrMsgTxt("windows size should be an scalar");
        }

        winSize = (int) mxGetScalar(in4);
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
	
	
	edrDist = mxGetPr(out1);
	edrCalculate(ts1, n1, ts2, n2, eps, winSize, dim1, edrDist);
    
    
    return;
    
}



