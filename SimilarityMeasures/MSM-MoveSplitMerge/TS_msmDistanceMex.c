#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// %     ************
// %
// %     Description
// %     ----------
// %     Compute the Move Split Merge between time series ts1 and ts2.
// % 
// %     Parameters
// %     ----------
// %     param ts1   :  n1 x dim, time series 1 matrix with the length of n1
// %     param ts2   :  n2 x dim, time series 2 matrix with the length of n2
// %
// %     Options
// %     ----------
// %     opt cost      : double, cost of Split/Merge operation.
// %     opt winSize   : integer, temporal constraint on the warping window
// %                              size. default value = -1
// %
// %     Returns
// %     -------
// %     msmDist      : double, The Move Split Merge distance between time series ts1 and ts2 
// %
// %     Other m-files required    : TS_euclideanDistance
// %     Subfunctions              : TS_msmCost
// %     MAT-files required        : none
// %     
// %     References
// %     ----------
// %     @article{stefan2013move,
// %       title={The move-split-merge metric for time series},
// %       author={Stefan, Alexandra and Athitsos, Vassilis and Das, Gautam},
// %       journal={IEEE transactions on Knowledge and Data Engineering},
// %       volume={25},
// %       number={6},
// %       pages={1425--1438},
// %       year={2013},
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

double msmCost(double *tsA, int nA, double *tsB, int nB, double *tsC, int nC, int dim, double cost, int idxA, int idxB, int idxC)
{
	double ax, ay, bx, by, cx, cy;
	double output, tmp1, tmp2;
	
	ax = tsA[idxA + nA * 0];
	ay = tsA[idxA + nA * 1];
	bx = tsB[idxB + nB * 0];
	by = tsB[idxB + nB * 1];
	cx = tsC[idxC + nC * 0];
	cy = tsC[idxC + nC * 1];
	
	if ( (bx <= ax && by <= ay && ax <= cx && ay <= cy) || (cx <= ax && cy <= ay && ax <= bx && ay <= by) )
	{
		output = cost;
	}
	else
	{
		tmp1 = euclideanCalculate(tsA, nA, tsB, nB, dim, idxA, idxB);
		tmp2 = euclideanCalculate(tsA, nA, tsC, nC, dim, idxA, idxC);
		
		if (tmp1 < tmp2)
		{
			output = cost + tmp1;
		}
		else
		{
			output = cost + tmp2;
		}
	}
	return output;
	
}

void msmCalculate(double *ts1, int n1, double *ts2, int n2, double cost, int winSize, int dim, double *msmDist)
{
	double D[n1][n2];
	int i, j;
	int jS, jF;
    int lenDiff;
	double tmp, tmp1, tmpMin;
	
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
	
	D[0][0] = euclideanCalculate(ts1, n1, ts2, n2, dim, 0, 0);
	
	for (i = 1; i < n1; i++)
	{
		D[i][0] = D[i-1][0] + msmCost(ts1, n1, ts1, n1, ts2, n2, dim, cost, i, i-1, 0);
	}
	
	for (j = 1; j < n2; j++)
	{
		D[0][j] = D[0][j-1] + msmCost(ts2, n2, ts1, n1, ts2, n2, dim, cost, j, 0, j-1);
	}
	
	for (i = 1; i < n1; i++)
	{
		if (winSize == -1)
		{
			jS = 1;
			jF = n2;
		}
		else
		{
			jS = i - winSize > 1 ? i - winSize : 1;
			jF = i + winSize < n2 ? i + winSize + 1 : n2;
		}
		
		for (j = jS; j < jF; j++)
		{
			tmp = euclideanCalculate(ts1, n1, ts2, n2, dim, i, j);
			D[i][j] = D[i-1][j-1] + tmp;
			
			tmp = D[i-1][j] + msmCost(ts1, n1, ts1, n1, ts2, n2, dim, cost, i, i-1, j);
			
			if (tmp < D[i][j])
			{
				D[i][j] = tmp;
			}
			
			tmp = D[i][j-1] + msmCost(ts2, n2, ts1, n1, ts2, n2, dim, cost, j, i, j-1);
			
			if (tmp < D[i][j])
			{
				D[i][j] = tmp;
			}
			

		}
	}
	
	
    msmDist[0] = D[n1-1][n2-1];
	
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

	
	/* Check the number of arguments */
	if(nrhs < 2 && nrhs > 4)
	{
		mexErrMsgTxt("Wrong number of input arguments.");
	}
	else if(nlhs != 1)
	{
		mexErrMsgTxt("Too many output arguments.");
	}
		
	/* Define local variables */
	double *ts1, *ts2;
	double cost;
	int winSize;
	int n1, dim1, n2, dim2;
	double *msmDist;
	
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
		cost = 0.1;
        winSize = -1;
    }
	else if(nrhs == 3)
    {
		cost = mxGetScalar(in3);
        winSize = -1;

    }
	else if(nrhs == 4)
    {
		cost = mxGetScalar(in3);
        winSize = (int) mxGetScalar(in4);
	}

	
	out1 = mxCreateDoubleMatrix( 1, 1, mxREAL);
	
	msmDist = mxGetPr(out1);
	
	msmCalculate(ts1, n1, ts2, n2, cost, winSize, dim1, msmDist);
    
    
    return;
    
}



