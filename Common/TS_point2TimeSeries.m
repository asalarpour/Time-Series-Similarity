function pt2trajDist = TS_point2TimeSeries(pt, ts)

%     ************
%
%     Description
%     ----------
%     Compute the shortest distance between queryPoint "pt" the timeseries "ts".
%     The Point to time series distance is the minimum of point-to-segment
%     distance between pt and all segment of ts
% 
%     Parameters
%     ----------
%     param pt      : 1 x dim, query point 
%     param ts    : m x dim, timeseries matrix with the length of m 
%
%     Returns
%     -------
%     pt2trajDist   :  double, Point-to-timeseries distance between pt and ts.
%
%     Other m-files required    : TS_point2SegmentDistance
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


totalSegDist = -ones(size(ts, 1) - 1, 1);

for i = 1: size(ts, 1) - 1
    totalSegDist(i) = TS_point2SegmentDistance(pt, ts(i, :) , ts(i + 1, :));
end

pt2trajDist = min(totalSegDist);

end
        