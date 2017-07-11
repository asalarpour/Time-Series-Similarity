function frechDist = TS_frechetDistance(ts1, ts2)

%     ************
%
%     Description
%     ----------
%     Compute the frechet distance from time series ts1 to time series ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Returns
%     -------
%     frechDist     :  double, the frechet distance between time series ts1 and ts2
%
%     Other m-files required    : TS_computeCriticalValues
%                                 TS_freeSpace
%                                 TS_reachableSpace
%     Subfunctions              : none
%     MAT-files required        : none
%     
%     References
%     ----------
%
%     Author
%     ----------
%     Amir Salarpour, Ph.D. Candidate, Artificial Intelligence
%     Bu-Ali Sina University, Hamedan, Iran, Dept. of Computer Engineering
%     email address : amir.salarpour@gmail.com  
%     Website       : http://www.salarpour.com
%     December 2016 : Last revision: 29-Jan-2017
%     
%     ************


critVal = TS_computeCriticalValues(ts1, ts2);
frechDist = critVal(1);

first = 1;
last = length(critVal);
middle = round( (first+ last) / 2);

while (first < last)
    eps = critVal(middle);
%     eps = 138;
    [freeSpace1, freeSpace2] = TS_freeSpace(ts1, ts2, eps);
    [rep, ~, ~] = TS_reachableSpace(freeSpace1, freeSpace2);
        
%     fprintf(' %d -----  %d ----- %f \n ', middle, rep, eps);
    
    if rep
        frechDist = eps;
        last = middle - 1;
    else
        first = middle + 1;
    end
    middle = round( (first+ last) / 2);
    
end
% frechDist = critVal;

end
