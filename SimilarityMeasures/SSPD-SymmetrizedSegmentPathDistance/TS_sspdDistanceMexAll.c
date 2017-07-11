#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>



double point2SegmentDistance(double *ts1, int n1, double *ts2, int n2, int dim, int idx1, int idx2)
{
	
	double A, B, C, D, E, F, G, H;
	double dotProd, lenSq, tmp, param, distStart, distEnd;
	double xProj, yProj;
	
	A = ts1[idx1 + n1 * 0] - ts2[idx2 + n2 * 0];
	B = ts1[idx1 + n1 * 1] - ts2[idx2 + n2 * 1];
	C = ts2[(idx2 + 1) + n2 * 0] - ts2[idx2 + n2 * 0];
	D = ts2[(idx2 + 1) + n2 * 1] - ts2[idx2 + n2 * 1];
    
	E = ts1[idx1 + n1 * 0] - ts2[(idx2 + 1) + n2 * 0];
	F = ts1[idx1 + n1 * 1] - ts2[(idx2 + 1) + n2 * 1];
	
	dotProd = A * C + B * D;
	lenSq = C * C + D * D;
	
	if (lenSq == 0)
	{
		tmp = sqrt(A * A + B * B);
		return tmp;
	}
	else
	{
		param = dotProd / lenSq;
		if (param < 0.00001 || param > 1)
		{
			
			distStart = sqrt(A * A + B * B);
			distEnd = sqrt(E * E + F * F);

            
			if (distStart > distEnd)
			{
				return distEnd;
			}
			else
			{
				return distStart;
			}
		}
		xProj = ts2[idx2 + n2 * 0] + param * C;
		yProj = ts2[idx2 + n2 * 1] + param * D;
        
		G = ts1[idx1 + n1 * 0] - xProj;
		H = ts1[idx1 + n1 * 1] - yProj;
		tmp = sqrt(G * G + H * H);
		return tmp;
	}
	
}

double sspdCalculate(double *ts1, int n1, double *ts2, int n2, int dim)
{
	int i, j;
	double tmp, tmpMin, tmpSum = 0;
	
	for (i = 0; i < n1; i++)
	{
		tmpMin = INFINITY;
		for (j = 0; j < n2 - 1; j++)
		{
			tmp = point2SegmentDistance(ts1, n1, ts2, n2, dim, i, j);
			if (tmp < tmpMin)
			{
				tmpMin = tmp;
			}
		}
		tmpSum += tmpMin;
	}
	
	tmpSum /= n1;
	
	return tmpSum;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Macros for the ouput and input arguments */
	#define out1 		plhs[0]
	#define in1 		prhs[0]

	
	/* Check the number of arguments */
	if(nrhs != 1)
	{
		mexErrMsgTxt("Wrong number of input arguments.");
	}
	else if(nlhs != 1)
	{
		mexErrMsgTxt("Bad output arguments.");
	}
		
    /* Define local variables */
	double 	*ts1, *ts2;
	int 	n1, dim1, n2, dim2;
	int		nFields, frstField = 0;
	int     NStructElems;
	int		i, j;
	mxArray    *tmp1, *tmp2;
	double *sspdDist;  
    double spd1, spd2;
    
	/* get input arguments */
    nFields = mxGetNumberOfFields(in1);
    NStructElems = mxGetNumberOfElements(in1);
	
	out1 = mxCreateDoubleMatrix( NStructElems, NStructElems, mxREAL);
	
	sspdDist = mxGetPr(out1);
	
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
			
            if (dim1 != dim2)
				mexErrMsgTxt("Two time series dimension must be the same");
            
            spd1 = sspdCalculate(ts1, n1, ts2, n2, dim1);
            spd2 = sspdCalculate(ts2, n2, ts1, n1, dim1);
    
            sspdDist[j + NStructElems * i] = (spd1 + spd2) / 2;
            
		}
	}
	
	mexPrintf("\n");	
	
    return;
    
}



