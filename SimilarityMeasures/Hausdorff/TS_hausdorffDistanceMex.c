#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// %     ************
// %
// %     Description
// %     ----------
// %     Compute the Hausdorff distance between time series ts1 and ts2.
// % 
// %     Parameters
// %     ----------
// %     param ts1   :  n1 x dim, time series 1 matrix with the length of n1
// %     param ts2   :  n2 x dim, time series 2 matrix with the length of n2
// %
// %
// %     Returns
// %     -------
// %     ehausDist       :  double, The Hausdorff distance between time series ts1 and ts2.
// %
// %     Other m-files required    : TS_directedHausdorff
// %     Subfunctions              : none
// %     MAT-files required        : none
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
// %     December 2016 : Last revision: 29-Jan-2017
// %     
// %     ************

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

double hausCalculate(double *ts1, int n1, double *ts2, int n2, int dim)
{
	int i, j;
	double tmp, tmpMin, tmpMax = 0;
	
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
		if ( tmpMin > tmpMax)
		{
			tmpMax = tmpMin;
		}
	}
	
	return tmpMax;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* Macros for the ouput and input arguments */
	#define out1 		plhs[0]
	#define in1 		prhs[0]
	#define in2 		prhs[1]

	
	/* Check the number of arguments */
	if(nrhs != 2)
	{
		mexErrMsgTxt("Wrong number of input arguments.");
	}
	else if(nlhs != 1)
	{
		mexErrMsgTxt("Bad output arguments.");
	}
		
	/* Define local variables */
	double *ts1, *ts2;
	int n1, dim1, n2, dim2;
	double *hausDist;
	double tmp1, tmp2;
	
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
	
	
	hausDist = mxGetPr(out1);
	
	tmp1 = hausCalculate(ts1, n1, ts2, n2, dim1);
    tmp2 = hausCalculate(ts2, n2, ts1, n1, dim1);
    
	if (tmp1 > tmp2)
	{
		hausDist[0] = tmp1;
	}
	else
	{
		hausDist[0] = tmp2;
	}
	
    return;
    
}



