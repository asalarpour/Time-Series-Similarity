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
	C = ts2[idx2 + 1 + n2 * 0] - ts2[idx2 + n2 * 0];
	D = ts2[idx2 + 1 + n2 * 1] - ts2[idx2 + n2 * 1];
    
	E = ts1[idx1 + n1 * 0] - ts2[(idx2 + 1) + n2 * 0];
	F = ts1[idx1 + n1 * 1] - ts2[(idx2 + 1) + n2 * 1];
	
	dotProd = A * C + B * D;
	lenSq = C * C + D * D;
	
	if (lenSq == 0)
	{
		tmp = sqrt(A * A + B * B);
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
				tmp = distEnd;
			}
			else
			{
				tmp = distStart;
			}
		}
		else
		{
			xProj = ts2[idx2 + n2 * 0] + param * C;
			yProj = ts2[idx2 + n2 * 1] + param * D;
			
			G = ts1[idx1 + n1 * 0] - xProj;
			H = ts1[idx1 + n1 * 1] - yProj;
			tmp = sqrt(G * G + H * H);
		}
		
	}
	//mexPrintf("%f  ", tmp);
	return tmp;
	
}

void criticValCalculate(double *ts1, int n1, double *ts2, int n2, int dim, double *criticVal, int *criticSize)
{
	int i, j;
	double tmp, tmp1, tmp2;
	double excepPoint;
	int position, count = 0;
	double Lij, Bij;
	
	tmp1 = 0;
	tmp2 = 0;
	for (i = 0; i < dim; i++)
	{
		tmp = ts1[0 + n1 * i] - ts2[0 + n2 * i];
		tmp1 += tmp * tmp;
		tmp = ts1[n1 - 1 + n1 * i] - ts2[n2 - 1 + n2 * i];
		tmp2 += tmp * tmp;
	}
	
	tmp1 = sqrt(tmp1);
	tmp2 = sqrt(tmp2);
	if (tmp1 > tmp2)
	{
		excepPoint = tmp1;
	}
	else
	{
		excepPoint = tmp2;
	}

	criticVal[count] = excepPoint;
	
    
	for (i = 0; i < n1 - 1; i++)
	{
		for (j = 0; j < n2 - 1; j++)
		{

			Lij = point2SegmentDistance(ts1, n1, ts2, n2, dim, i, j);
			if (Lij > excepPoint)
			{
				count++;
				criticVal[count] = Lij;
			}
			Bij = point2SegmentDistance(ts2, n2, ts1, n1, dim, j, i);
			if (Bij > excepPoint)
			{
				count++;
				criticVal[count] = Bij;
			}

		}
	}
    
	
	*criticSize = count;


	i = 0;
	while ( i <= *criticSize)
	{
        position = i;
		j = i + 1;
		while (j <= *criticSize)
		{
//			if ( fabs(criticVal[i] - criticVal[j]) <= 0.0017 )
			if ( criticVal[i] == criticVal[j] )
			{
				criticVal[j] = criticVal[*criticSize];
				*criticSize = *criticSize - 1;
			}
			if (criticVal[position] > criticVal[j])
			{
				position = j;
			}
			
			j++;
		}
        tmp = criticVal[i];
        criticVal[i] = criticVal[position];
        criticVal[position] = tmp;
		i++;
	}
	
	

	
}

void circleLineIntersectCalculation(double *ts1, int n1, double *ts2, int n2, int dim, double radius, int idx1, int idx2, double intersect[][dim])
{
	double xCenter, yCenter;
	double lp1x, lp1y, lp2x, lp2y;
	double rac;
	double x, y, x1, x2, y1, y2;
	double m, b, A, B, C, delta, sDelta;
	
	xCenter = ts1[idx1 + n1 * 0];
	yCenter = ts1[idx1 + n1 * 1];
	lp1x = ts2[idx2 + n2 * 0];
	lp1y = ts2[idx2 + n2 * 1];
	lp2x = ts2[idx2 + 1 + n2 * 0];
	lp2y = ts2[idx2 + 1 + n2 * 1]; 
	
	if ( lp1x == lp2x)
	{
		rac = sqrt (radius * radius) - (lp1x - xCenter) * (lp1x - xCenter);
		y1 = yCenter + rac;
		y2 = yCenter - rac;
		intersect[0][0] = lp1x;
		intersect[0][1] = y1;
		intersect[1][0] = lp1x;
		intersect[1][1] = y2;
	}
	else
	{
		m = (lp2y - lp1y) / (lp2x - lp1x);
		b = lp2y - m * lp2x;
		A = m * m + 1;
		B = 2 * (m * b - m * yCenter - xCenter);
		C = (yCenter * yCenter) - (radius * radius) + (xCenter * xCenter) - (2 * b * yCenter) + (b * b);
		delta = B * B - 4 * A * C;
		
		if (delta <= 0)
		{
			x = -B / (2 * A);
			y = m * x + b;
			intersect[0][0] = x;
			intersect[0][1] = y;
			intersect[1][0] = x;
			intersect[1][1] = y;
		}
		else
		{
			if (delta > 0)
			{
				sDelta = sqrt(delta);
				x1 = (-B + sDelta) / (2 * A);
				y1 = m * x1 + b;
				x2 = (-B - sDelta) / (2 * A);
				y2 = m * x2 + b;
				intersect[0][0] = x1;
				intersect[0][1] = y1;
				intersect[1][0] = x2;
				intersect[1][1] = y2;
			}
			else
			{
				mexErrMsgTxt("The intersection between circle and line is supposed to exist");
			}
		}
	}
}

void freelineCalculate(double *ts1, int n1, double *ts2, int n2, int dim, double eps, int idx1, int idx2, double *valPt1, double *valPt2)
{
	double px, py;
	double seg1x, seg1y, seg2x, seg2y;
	double tmp, tmp1, tmp2;
	double segLen, segLenS;
	double intersect[2][dim];
	double i1x, i1y, u1, i2x, i2y, u2;
	double tmpArray[4];
	int position;
	int i, j;
	double val1, val2;
	
 	px = ts1[idx1 + n1 * 0];
	py = ts1[idx1 + n1 * 1];
	seg1x = ts2[idx2 + n2 * 0];
	seg1y = ts2[idx2 + n2 * 1];
	seg2x = ts2[idx2 + 1 + n2 * 0];
	seg2y = ts2[idx2 + 1 + n2 * 1]; 
	


	if (seg1x == seg2x && seg1y == seg2y)
	{
 		tmp1 = px - seg1x;
		tmp1 = py - seg1y;
		tmp = sqrt(tmp1 * tmp1 + tmp2 * tmp2);
		if (tmp > eps)
		{
			val1 = -1;
			val2 = -1;
		}
		else
		{
			val1 = 0;
			val2 = 1;
		} 
	}
	else
	{
		tmp = point2SegmentDistance(ts1, n1, ts2, n2, dim, idx1, idx2);
		//mexPrintf(" %f ", tmp);
		//mexPrintf(" %f+%f+%f+%f+%f+%f+%f  ", tmp, px, py, seg1x, seg1y, seg2x, seg2y);
//		mexPrintf("%d +++ %d ", idx1, idx2);

		if (tmp > eps)
		{
 			val1 = -1;
			val2 = -1; 
		}
		else
		{
 			tmp1 = seg1x - seg2x;
			tmp2 = seg1y - seg2y;
			segLen = sqrt(tmp1 * tmp1 + tmp2 * tmp2);
			segLenS = segLen * segLen; 
 			circleLineIntersectCalculation(ts1, n1, ts2, n2, dim, eps, idx1, idx2, intersect);
			

 			if (intersect[0][0] != intersect[1][0] || intersect[0][1] != intersect[1][1])
			{
				i1x = intersect[0][0];
				i1y = intersect[0][1];
				u1 = ((i1x - seg1x) * (seg2x - seg1x) + (i1y - seg1y) * (seg2y - seg1y)) / segLenS;
				
				i2x = intersect[1][0];
				i2y = intersect[1][1];
				u2 = ((i2x - seg1x) * (seg2x - seg1x) + (i2y - seg1y) * (seg2y - seg1y)) / segLenS;
				
				tmpArray[0] = 0;
				tmpArray[1] = 1;
				tmpArray[2] = u1;
				tmpArray[3] = u2;
				
				for (i = 0; i < 4; i++)
				{
					position = i;
					for (j = i + 1; j < 4; j++)
					{
						if (tmpArray[position] > tmpArray[j] )
						{
							position = j;
						}
					}
					tmp = tmpArray[i];
					tmpArray[i] = tmpArray[position];
					tmpArray[position] = tmp;
				}
				
				val1 = tmpArray[1];
				val2 = tmpArray[2];
			}
			else
			{
				if ( px == seg1x && py == seg1y)
				{
					val1 = 0;
					val2 = 0;
				}
				else 
				{
					if(px == seg2x && py == seg2y)
					{
						val1 = 1;
						val2 = 1;
					}
					else
					{
						i1x = intersect[0][0];
						i1y = intersect[0][1];
						u1 = ((i1x - seg1x) * (seg2x - seg1x) + (i1y - seg1y) * (seg2y - seg1y)) / segLenS;
						
						if ( u1 >= 0 && u1 <= 1)
						{
							val1 = u1;
							val2 = u1;
						}
						else
						{
							val1 = -1;
							val2 = -1;
						}
					}
				}
				
			} 
			
		}
	}
	*valPt1 = val1;
	*valPt2 = val2;
}

void freespaceCalculate(double *ts1, int n1, double *ts2, int n2, int dim, double eps, double LF[][n2][dim], double BF[][n2][dim])
{
	double tmp1, tmp2;
	int i, j;
	for (j = 0; j < n2; j++)
	{
		for (i = 0; i < n1 - 1; i++)
		{
			freelineCalculate(ts2, n2, ts1, n1, dim, eps, j, i, &tmp1, &tmp2);
			LF[i][j][0] = tmp1;
			LF[i][j][1] = tmp2;
			
		}
		
	}
	
	for (j = 0; j < n2 - 1; j++)
	{
		for (i = 0; i < n1; i++ )
		{
			freelineCalculate(ts1, n1, ts2, n2, dim, eps, i, j, &tmp1, &tmp2);
			BF[i][j][0] = tmp1;
			BF[i][j][1] = tmp2;
			//mexPrintf("[%f,%f] ", tmp1, tmp2);
		}
		//mexPrintf("\n");
	}
	
}

bool reachableCalculate(double *ts1, int n1, double *ts2, int n2, int dim, double eps)
{
	
	double LF[n1][n2][dim], BF[n1][n2][dim];
	bool LR[n1][n2], BR[n1][n2];
	bool rep;
	int i, j;
	
	freespaceCalculate(ts1, n1, ts2, n2, dim, eps, LF, BF);
	
	if (!(LF[0][0][0] == 0 && BF[0][0][0] == 0 &&  LF[n1-2][n2-1][1] == 1 && BF[n1-1][n2-2][1] == 1 ))
	{
		rep = false;
		return rep;
	}
	else
	{
		LR[0][0] = true;
		BR[0][0] = true;

		for (i = 1; i < n1 - 1; i++)
		{
			if (LF[i][0][0] != -1 && LF[i][0][1] != -1 && LF[i-1][0][0] == 0 && LF[i-1][0][1] == 1 )
			{
				LR[i][0] = true;
			}
			else
			{
				LR[i][0] = false;
			}

		}
		for (j = 1; j < n2 - 1; j++)
		{
			if (BF[0][j][0] != -1 && BF[0][j][1] != -1 && BF[0][j-1][0] == 0 && BF[0][j-1][1] == 1 )
			{
				//mexPrintf("%d    %f,%f   %f,%f \n",j,  BF[0][j][0], BF[0][j][1],  BF[0][j-1][0], BF[0][j-1][1]);
				BR[0][j] = true;
			}
			else
			{
				BR[0][j] = false;
			}

		}
		
		for (i = 0; i < n1 - 1; i++)
		{
			for (j = 0; j < n2 - 1; j++)
			{
				
				if (LR[i][j] || BR[i][j])
				{
					if (LF[i][j+1][0] != -1 && LF[i][j+1][1] != -1) 
					{
						LR[i][j+1] = true;
					}
					else
					{
						LR[i][j+1] = false;
					}
					
					if (BF[i+1][j][0] != -1 && BF[i+1][j][1] != -1) 
					{
						BR[i+1][j] = true;
					}
					else
					{
						BR[i+1][j] = false;
					}	
					
				}
				else
				{
					LR[i][j+1] = false;
					BR[i+1][j] = false;
				}
			}
		}
/* 		for (i = 0; i< n1-1; i++)
		{
			mexPrintf("%d ", i);
			for (j = 0; j < n2 - 1; j++)
			{
				mexPrintf("%d,%d  ", BR[i][j], LR[i][j]);
			}
			mexPrintf("\n");
		} */
		
/* 		mexPrintf("%d,%d  ", BR[n1-2][n2-2], LR[n1-2][n2-2]);
 */		
		rep = BR[n1-2][n2-2] || LR[n1-2][n2-2];		
	}
	return rep;
}

double frechetCalculate(double *ts1, int n1, double *ts2, int n2, int dim)
{
	double criticVal[n1 * n2 * 2];

	int criticSize = 0;
	double eps, outDist;
	int first, last, middle;
	bool reach;
	
	criticValCalculate(ts1, n1, ts2, n2, dim, criticVal, &criticSize);
	
	first = 0;
	last = criticSize;
	middle = (first + last) / 2;
	outDist = criticVal[0];
	

	while (first <= last)
	{
		eps = criticVal[middle];
		//eps = 138;
		
		reach = reachableCalculate(ts1, n1, ts2, n2, dim, eps);
		
		//mexPrintf("%d -----  %d ----- %f \n ", middle, reach, eps);
		if (reach)
		{
			last = middle - 1;
			outDist = eps;
			
		}
		else
		{
			first = middle + 1;
		}
		middle = (first + last)/2;
		
	}
	
	return outDist;
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
	double *frechetDist;
    
	/* get input arguments */
    nFields = mxGetNumberOfFields(in1);
    NStructElems = mxGetNumberOfElements(in1);
	
	out1 = mxCreateDoubleMatrix( NStructElems, NStructElems, mxREAL);
	
	frechetDist = mxGetPr(out1);    
	
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
            
            frechetDist[j + NStructElems * i] = frechetCalculate(ts1, n1, ts2, n2, dim1);

		}
	}
	
	mexPrintf("\n");	
    
    return;
    
}



