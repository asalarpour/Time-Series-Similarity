function hausDist = TS_hausdorffDistance(ts1, ts2)

%     ************
%
%     Description
%     ----------
%     Compute the Hausdorff distance between time series ts1 and ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%
%     Returns
%     -------
%     ehausDist       :  double, The Hausdorff distance between time series ts1 and ts2.
%
%     Other m-files required    : TS_directedHausdorff
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




hausDist = max([TS_directedHausdorff(ts1,ts2),...
    TS_directedHausdorff(ts2,ts1)]);

end