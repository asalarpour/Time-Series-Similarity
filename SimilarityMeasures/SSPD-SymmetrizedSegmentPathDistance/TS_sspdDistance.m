function sspdDist = TS_sspdDistance(ts1, ts2)

%     ************
%
%     Description
%     ----------
%     Compute the sspd-distance between time series ts1 and ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Returns
%     -------
%     sspdDist      :  double, sspd-distance between time series ts1 and ts2
%
%     Other m-files required    : TS_spdDistance
%     Subfunctions              : none
%     MAT-files required        : none
%     
%     References
%     ----------
%     @techreport{eiter1994computing,
%       title={Computing discrete Fr{\'e}chet distance},
%       author={Eiter, Thomas and Mannila, Heikki},
%       year={1994},
%       institution={Citeseer}
%     }
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


sspdDist = TS_spdDistance(ts1, ts2) + TS_spdDistance(ts2, ts1);
sspdDist = sspdDist/ 2;

end