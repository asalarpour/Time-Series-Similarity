#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// %     ************
// %
// %     Description
// %     ----------
// %     Compute the sspd-distance between time series ts1 and ts2.
// % 
// %     Parameters
// %     ----------
// %     param ts1   :  m x dim, time series 1 matrix with the length of m
// %     param ts2   :  n x dim, time series 2 matrix with the length of n
// %
// %     Returns
// %     -------
// %     sspdDist      :  double, sspd-distance between time series ts1 and ts2
// %
// %     Other m-files required    : TS_spdDistance
// %     Subfunctions              : none
// %     MAT-files required        : none
// %     
// %     References
// %     ----------
// %     @techreport{eiter1994computing,
// %       title={Computing discrete Fr{\'e}chet distance},
// %       author={Eiter, Thomas and Mannila, Heikki},
// %       year={1994},
// %       institution={Citeseer}
// %     }
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
	double *sspdDist;
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
	
	
	sspdDist = mxGetPr(out1);
	
	tmp1 = sspdCalculate(ts1, n1, ts2, n2, dim1);
    tmp2 = sspdCalculate(ts2, n2, ts1, n1, dim1);
    
	sspdDist[0] = (tmp1 + tmp2) / 2;
	
    return;
    
}



