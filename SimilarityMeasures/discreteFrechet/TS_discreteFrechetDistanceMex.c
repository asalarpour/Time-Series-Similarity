#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// %     ************
// %
// %     Description
// %     ----------
// %     Compute the discret frechet distance between time series ts1 and ts2.
// % 
// %     Parameters
// %     ----------
// %     param ts1   :  n1 x dim, time series 1 matrix with the length of n1
// %     param ts2   :  n2 x dim, time series 2 matrix with the length of n2
// %
// %     Returns
// %     -------
// %     disFrechDist  :  double, the discret frechet distance between time series ts1 and ts2
// %
// %     Other m-files required    : TS_euclideanDistance
// %     Subfunctions              : TS_func
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

double recursiveFunc(double *ts1, int n1, double *ts2, int n2, int dim, double *ca, int idx1 , int idx2)
{
	double output;
	double tmp, tmp1, tmp2;
	int i;
	
	if (ca[idx1 + n1 * idx2] > -1)
	{
		output = ca[idx1 + n1 * idx2];
		return output;
	}
	else if ( idx1 == 0 && idx2 == 0)
	{
		tmp1 = 0;
		for ( i = 0; i < dim; i++)
		{
			tmp2 = ts1[idx1 + n1 * i] - ts2[idx2 + n2 * i];
			tmp1 +=  tmp2 * tmp2;
		}
		
		ca[idx1 + n1 * idx2] = sqrt(tmp1);
	}
	else if ( idx1 > 0 && idx2 == 0)
	{
		tmp = recursiveFunc(ts1, n1, ts2, n2, dim, ca, idx1-1 , idx2);
		tmp1 = 0;
		for ( i = 0; i < dim; i++)
		{
			tmp2 = ts1[idx1 + n1 * i] - ts2[idx2 + n2 * i];
			tmp1 +=  tmp2 * tmp2;
		}
		if (tmp > sqrt(tmp1))
		{
			ca[idx1 + n1 * idx2] = tmp;
		}
		else
		{
			ca[idx1 + n1 * idx2] = sqrt(tmp1);
		}
	}
	else if ( idx1 == 0 && idx2 > 0)
	{
		tmp = recursiveFunc(ts1, n1, ts2, n2, dim, ca, idx1 , idx2-1);
		tmp1 = 0;
		for ( i = 0; i < dim; i++)
		{
			tmp2 = ts1[idx1 + n1 * i] - ts2[idx2 + n2 * i];
			tmp1 +=  tmp2 * tmp2;
		}
		if (tmp > sqrt(tmp1))
		{
			ca[idx1 + n1 * idx2] = tmp;
		}
		else
		{
			ca[idx1 + n1 * idx2] = sqrt(tmp1);
		}
	}
		else if ( idx1 > 0 && idx2 > 0)
	{
		tmp = recursiveFunc(ts1, n1, ts2, n2, dim, ca, idx1-1, idx2-1);
		tmp1 = recursiveFunc(ts1, n1, ts2, n2, dim, ca, idx1-1, idx2);
		tmp2 = recursiveFunc(ts1, n1, ts2, n2, dim, ca, idx1, idx2-1);
		if ( tmp1 < tmp)
		{
			tmp = tmp1;
		}
		if (tmp2 < tmp)
		{
			tmp = tmp2;
		}
		
		tmp1 = 0;
		for ( i = 0; i < dim; i++)
		{
			tmp2 = ts1[idx1 + n1 * i] - ts2[idx2 + n2 * i];
			tmp1 +=  tmp2 * tmp2;
		}
		if (tmp > sqrt(tmp1))
		{
			ca[idx1 + n1 * idx2] = tmp;
		}
		else
		{
			ca[idx1 + n1 * idx2] = sqrt(tmp1);
		}
	}
	output = ca[idx1 + n1 * idx2];
	return output;
}


double disFrechetCalculate(double *ts1, int n1, double *ts2, int n2, int dim)
{
	double ca[n1 * n2];
	int i, j;
	double tmpDist;
	
	for ( i = 0; i < n1; i++)
	{
		for ( j = 0; j < n2; j++)
		{
			ca[i + n1 * j] = -1;
		}
	}
	
	
	tmpDist = recursiveFunc(ts1, n1, ts2, n2, dim, ca, n1 -1 , n2 -1);
	
	return tmpDist;
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
	double *disFrechetDist;
	
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
	
	
	disFrechetDist = mxGetPr(out1);
	
	disFrechetDist[0] = disFrechetCalculate(ts1, n1, ts2, n2, dim1);
	
    return;
    
}



