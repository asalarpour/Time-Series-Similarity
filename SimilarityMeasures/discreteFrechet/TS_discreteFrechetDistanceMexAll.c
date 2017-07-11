#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


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
	double *disFrechetDist;
	
	/* get input arguments */
    nFields = mxGetNumberOfFields(in1);
    NStructElems = mxGetNumberOfElements(in1);
	
	out1 = mxCreateDoubleMatrix( NStructElems, NStructElems, mxREAL);
	
	disFrechetDist = mxGetPr(out1);
	
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
            
            disFrechetDist[j + NStructElems * i] = disFrechetCalculate(ts1, n1, ts2, n2, dim1);

		}
	}
	
    return;
    
}



