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

void lcssCalculate(double *ts1, int n1, double *ts2, int n2, double eps, int winSize, int dim, double *lcssDist)
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
		tmp = (double) n2;
	}
	else
	{
		int lenDiff = n2 - n1;
		tmp = (double) n1;
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
				D[i][j] = D[i-1][j-1] + 1;
			}
			else
			{
				if (D[i][j-1] > D[i-1][j])
				{
					D[i][j] = D[i][j-1];
				}
				else
				{
					D[i][j] = D[i-1][j];
				}
			}
			
		}
	}
	
    
    lcssDist[0] = 1 - (D[n1][n2] / tmp);
	
	return;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Macros for the ouput and input arguments */
	#define out1 		plhs[0]
	#define in1 		prhs[0]
	#define in2 		prhs[1]
	#define in3 		prhs[2]

	
	/* Check the number of arguments */
	if(nrhs != 1 && nrhs != 2 && nrhs != 3)
	{
		mexErrMsgTxt("Wrong number of input arguments.");
	}
	else if(nlhs != 1)
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
	double 	*lcssDist;
	double eps;	
	
	/* check to make sure winSize is a scalar */
	if(nrhs == 1)
    {
		eps = INFINITY;
        winSize = -1;
    }
	else if(nrhs == 2)
    {
        eps = mxGetScalar(in2);
		winSize = -1;
    }
	else if(nrhs == 3)
    {
        eps = mxGetScalar(in2);

        winRatio = (int) mxGetScalar(in3);
    }
	
	/* get input arguments */
    nFields = mxGetNumberOfFields(in1);
    NStructElems = mxGetNumberOfElements(in1);
	
	out1 = mxCreateDoubleMatrix( NStructElems, NStructElems, mxREAL);
	
	lcssDist = mxGetPr(out1);
	
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
				if (winSize != -1)
					winSize = (int) tmp + n2 * winRatio;
			}
			else
			{
				if (winSize != -1)
					winSize = (int) tmp + n1 * winRatio;
			}
			
			if (dim1 != dim2)
				mexErrMsgTxt("Two time series dimension must be the same");
			lcssCalculate(ts1, n1, ts2, n2, eps, winSize, dim1, &lcssDist[j + NStructElems * i]);
			

		}
	}
	

	mexPrintf("\n");
    
    
    return;
    
}



