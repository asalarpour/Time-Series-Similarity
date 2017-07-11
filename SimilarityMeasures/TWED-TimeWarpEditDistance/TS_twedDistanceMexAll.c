#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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
    #define out2 		plhs[1]
	#define in1 		prhs[0]
	#define in2 		prhs[1]
	#define in3 		prhs[2]
	#define in4 		prhs[3]

	
	/* Check the number of arguments */
	if(nrhs < 1 && nrhs > 4)
	{
		mexErrMsgTxt("Wrong number of input arguments.");
	}
	else if(nlhs != 1 && nlhs != 2)
	{
		mexErrMsgTxt("Too many output arguments.");
	}
		
    /* Define local variables */
	double 	winRatio, tmp;
	int 	winSize = -1;
	double 	*ts1, *ts2;
	int 	n1, dim1, n2, dim2;
	int		nFields, frstField = 0;
	int     NStructElems;
	int		i, j;
	mxArray    *tmp1, *tmp2;
	double 	*twedDist, *twedLenMax;
    double lambda, nu;
    
    /* check to make sure winSize is a scalar */
	if(nrhs == 1)
    {
		lambda = 1;
		nu = 0;
        winSize = -1;
    }
	else if(nrhs == 2)
    {
		lambda = mxGetScalar(in2);
		nu = 0;
        winSize = -1;
    }
	else if(nrhs == 3)
    {
        lambda = mxGetScalar(in2);
		nu = mxGetScalar(in3);
        winSize = -1;
    } 
	else if(nrhs == 4)
    {
        lambda = mxGetScalar(in2);
		nu = mxGetScalar(in3);
        winRatio = (int) mxGetScalar(in4);
    }    
    

    /* get input arguments */
    nFields = mxGetNumberOfFields(in1);
    NStructElems = mxGetNumberOfElements(in1);
	
	out1 = mxCreateDoubleMatrix( NStructElems, NStructElems, mxREAL);
	out2 = mxCreateDoubleMatrix( NStructElems, NStructElems, mxREAL);
	
	twedDist = mxGetPr(out1);
	twedLenMax = mxGetPr(out2);
    
    double percent;
	mexPrintf("\n");
	for(i = 0; i < NStructElems; i++)
	{
		percent = (double) 100 * (i + 1) / NStructElems;
		mexPrintf("Percent Done : %3.1f \n", percent);
		mexEvalString("drawnow");		

		
		for(j = i + 1; j < NStructElems; j++)
		{
			tmp1 = mxGetFieldByNumber(in1, i, frstField);
			tmp2 = mxGetFieldByNumber(in1, j, frstField);
			
			ts1 = mxGetPr(tmp1);
			ts2 = mxGetPr(tmp2);
			
			n1 = mxGetM(tmp1);
			dim1 = mxGetN(tmp1);
			
			n2 = mxGetM(tmp2);
			dim2 = mxGetN(tmp2);
			
			tmp = abs(n2 - n1);
			
			if (n1 > n2)
			{
				twedLenMax[j + NStructElems * i] = n1;
				if (winSize != -1)
					winSize = (int) tmp + n2 * winRatio;
			}
			else
			{
				twedLenMax[j + NStructElems * i] = n2;
				if (winSize != -1)
					winSize = (int) tmp + n1 * winRatio;
			}
			
			if (dim1 != dim2)
				mexErrMsgTxt("Two time series dimension must be the same");
			twedCalculate(ts1, n1, ts2, n2, lambda, nu, winSize, dim1, &twedDist[j + NStructElems * i]);
			

		}
	}
    
    
    return;
    
}



