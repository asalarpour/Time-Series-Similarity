function critVal = TS_computeCriticalValues(ts1, ts2)

%     ************
%
%     Description
%     ----------
%     Compute all the critical values between time series ts1 and ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Returns
%     -------
%     critVal       : 1 x k, all critical values between time series ts1 and ts1
%
%     Other m-files required    : TS_euclideanDistance, TS_point2SegmentDistance
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
n = size(ts2, 1);


distOrigin = TS_euclideanDistance(ts1(1, :), ts2(1, :));
distEnd = TS_euclideanDistance(ts1(end, :), ts2(end, :));

excepPoint = max([distOrigin, distEnd]);

critVal = excepPoint;

for i = 1: m - 1
    for j = 1: n - 1
        
        Lij = TS_point2SegmentDistance(ts2(j, :), ts1(i, :), ts1(i + 1, :));
        
        if Lij > excepPoint
            critVal = [critVal Lij]; %#ok<AGROW>
        end
        
        Bij = TS_point2SegmentDistance(ts1(i, :), ts2(j, :), ts2(j + 1, :));
        
        if Bij > excepPoint
            critVal = [critVal Bij]; %#ok<AGROW>
        end
        
    end
end
critVal = sort(unique(critVal));

end