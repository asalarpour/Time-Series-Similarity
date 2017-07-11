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
    #define out2 		plhs[1]
	#define in1 		prhs[0]
	#define in2 		prhs[1]
	#define in3 		prhs[2]

	
	/* Check the number of arguments */
	if(nrhs < 1 && nrhs > 3)
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
	double 	*msmDist, *msmLenMax;
    double cost;

    /* check to make sure winSize is a scalar */
	if(nrhs == 1)
    {
        cost = 0.1;
        winSize = -1;
    }
	else if(nrhs == 2)
    {
        cost = mxGetScalar(in2);
		winSize = -1;
    }
	else if(nrhs == 3)
    {
        cost = mxGetScalar(in2);
        winRatio = (int) mxGetScalar(in3);
    }	
    
    /* get input arguments */
    nFields = mxGetNumberOfFields(in1);
    NStructElems = mxGetNumberOfElements(in1);
	
	out1 = mxCreateDoubleMatrix( NStructElems, NStructElems, mxREAL);
	out2 = mxCreateDoubleMatrix( NStructElems, NStructElems, mxREAL);
	
	msmDist = mxGetPr(out1);
	msmLenMax = mxGetPr(out2);
    
    
    
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
				msmLenMax[j + NStructElems * i] = n1;
				if (winSize != -1)
					winSize = (int) tmp + n2 * winRatio;
			}
			else
			{
				msmLenMax[j + NStructElems * i] = n2;
				if (winSize != -1)
					winSize = (int) tmp + n1 * winRatio;
			}
			
			if (dim1 != dim2)
				mexErrMsgTxt("Two time series dimension must be the same");
			msmCalculate(ts1, n1, ts2, n2, cost, winSize, dim1, &msmDist[j + NStructElems * i]);
			

		}
	}
    
    return;
    
}



