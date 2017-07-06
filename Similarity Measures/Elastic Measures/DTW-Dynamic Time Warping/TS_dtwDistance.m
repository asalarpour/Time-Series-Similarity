function [dtwDist, dtwLen] = TS_dtwDistance(ts1, ts2, winSize)
    
%     ************
%
%     Description
%     ----------
%     Compute the Dynamic-Time Warping distance between time series ts1 and ts2.
% 
%     Parameters
%     ----------
%     param ts1   :   m x dim, time series 1 matrix with the length of m
%     param ts2   :   n x dim, time series 2 matrix with the length of n
%
%     Options
%     ----------
%     opt winSize   :   integer,   temporal constraint on the warping window
%                              size. default value = -1 
%
%     Returns
%     -------
%     dtwDist       :   double, The Dynamic-Time Warping distance between time series ts1 and ts2.
%     dtwLen        :   integer, Show the length of warping path between time series
%
%     Other m-files required    : TS_euclideanDistance
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
%     December 2016 : Last revision: 27-Jan-2017
%     
%     ************

if ~exist('winSize','var')
    winSize = -1;
end

if size(ts1, 2) == size(ts2, 2)
    dim = size(ts2, 2); %#ok<NASGU>
else
    error('Two time series dimension must be the same')
end

m = size(ts1, 1);
n = size(ts2, 1);

D = Inf(m + 1, n + 1);
L = zeros(m + 1, n + 1);

D(1, 1) = 0;

if winSize ~= -1
    winSize = max([winSize, abs(n - m)]);
end


for i = 2: m + 1
    if winSize == -1
        jS = 2;
        jF = n + 1;
    else
        jS = max([2, i - winSize]);
        jF = min ([n + 1, i + winSize]);
    end
    for j = jS:jF
        
        tmpDist = TS_euclideanDistance(ts1(i - 1, :), ts2(j - 1, :));
        minVal = min([D(i, j - 1), D(i - 1, j - 1), D(i - 1, j)]);
        D(i, j) = tmpDist + minVal;
        L(i, j) = min([L(i, j - 1), L(i - 1, j - 1), L(i - 1, j)]) + 1;

    end
end


dtwDist = D(m + 1, n + 1);
dtwLen = L(m + 1, n + 1);

end
