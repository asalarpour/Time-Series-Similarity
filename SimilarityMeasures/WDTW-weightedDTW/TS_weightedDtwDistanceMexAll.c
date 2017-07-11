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
	#define out3 		plhs[2]    
	#define in1 		prhs[0]
	#define in2 		prhs[1]
	#define in3 		prhs[2]
	#define in4 		prhs[3]
	
	/* Check the number of arguments */
	if(nrhs < 2 && nrhs > 4)
	{
		mexErrMsgTxt("Wrong number of input arguments.");
	}
	else if(nlhs != 1 && nlhs != 2 && nlhs != 3)
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
	double 	*wdtwDist, *wdtwLen, *wdtwLenMax;
    double g, weightMax;

    
    
    /* check to make sure winSize is a scalar */
	if(nrhs == 1)
    {
		g = 0.25;
		weightMax = 1;
        winSize = -1;
    }
	else if(nrhs == 2)
    {
        g = mxGetScalar(in2);
		weightMax = 1;
        winSize = -1;
    }
    else if(nrhs == 3)
    {
        g = mxGetScalar(in2);
		weightMax = mxGetScalar(in3);
        winSize = -1;
    }
	else if(nrhs == 4)
    {
        g = mxGetScalar(in2);
		weightMax = mxGetScalar(in3);
        winRatio = mxGetScalar(in4);
    }
    
    /* get input arguments */
    nFields = mxGetNumberOfFields(in1);
    NStructElems = mxGetNumberOfElements(in1);
	
	out1 = mxCreateDoubleMatrix( NStructElems, NStructElems, mxREAL);
	out2 = mxCreateDoubleMatrix( NStructElems, NStructElems, mxREAL);
	out3 = mxCreateDoubleMatrix( NStructElems, NStructElems, mxREAL);
	
	wdtwDist = mxGetPr(out1);
	wdtwLen = mxGetPr(out2);
	wdtwLenMax = mxGetPr(out3);
	
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
				wdtwLenMax[j + NStructElems * i] = n1;
				if (winSize != -1)
					winSize = (int) tmp + n2 * winRatio;
			}
			else
			{
				wdtwLenMax[j + NStructElems * i] = n2;
				if (winSize != -1)
					winSize = (int) tmp + n1 * winRatio;
			}
			
			if (dim1 != dim2)
				mexErrMsgTxt("Two time series dimension must be the same");
			wdtwCalculate(ts1, n1, ts2, n2, g, weightMax, winSize, dim1, &wdtwDist[j + NStructElems * i] , &wdtwLen[j + NStructElems * i]);

		}
	}    
    
    return;
    
}



