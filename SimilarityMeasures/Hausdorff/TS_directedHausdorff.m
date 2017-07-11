function dirHausDist = TS_directedHausdorff(ts1, ts2)

%     ************
%
%     Description
%     ----------
%     Compute the directed hausdorff distance from time series ts1 to time series ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Returns
%     -------
%     dirHausDist   :  double, directed hausdorff from time series ts1 to time series ts2
%
%     Other m-files required    : TS_point2time series
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

    

m = size(ts1, 1);
n = size(ts2, 1); %#ok<NASGU>
  
tmpHaus = -ones(m, 1);

for i = 1: m
    tmpHaus(i) = TS_point2TimeSeries( ts1(i, :), ts2);
end
dirHausDist = max(tmpHaus);

end
