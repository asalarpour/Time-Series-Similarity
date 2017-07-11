function spdDist = TS_spdDistance(ts1, ts2)

%     ************
%
%     Description
%     ----------
%     Compute the spd-distance of time series ts2 from time series ts1
%     The spd-distance is the sum of the all the point-to-time series
%     distance of points of ts1 from time series ts2
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Returns
%     -------
%     spdDist       :  double, spd-distance of time series ts2 from time series ts1
%
%     Other m-files required    : TS_point2TimeSeries
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


m = size(ts1, 1);
n = size(ts2, 1); %#ok<NASGU>

tmpSpd = [];

for i = 1: m
    tmp = TS_point2TimeSeries(ts1(i, :), ts2);
    tmpSpd = [tmpSpd tmp]; %#ok<AGROW>
end

spdDist = sum(tmpSpd) / length(tmpSpd);

end
