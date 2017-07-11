function [wdtwDist, wdtwLen] = TS_weightedDtwDistance(ts1, ts2, g, weightMax, winSize)
    

%     ************
%
%     Description
%     ----------
%     Compute the Weighted Dynamic-Time Warping distance between time series ts1 and ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Options
%     ----------
%     opt winSize   :   integer,   temporal constraint on the warping window
%                              size. default value = -1
%     opt g         :   double,  empirical constant that controls the curvature (slope) of the function
%                              g = 0 it is constant, g = 0.05 it is linear, g = 0.25 it is sigmoid, g = 3 two distinct weight
%                              the range of optimal g is between 0.01 to 0.6
%     opt weightMax : double,  maximum weight to normalze
%
%     Returns
%     -------
%     wdtwDist       : double, The Weighted Dynamic-Time Warping distance between time series ts1 and ts2
%     wdtwLen        :   integer, Show the length of warping path between time series
%
%     Other m-files required    : TS_euclideanDistance
%     Subfunctions              : TS_weightFun
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

if ~exist('g','var')
    g = 0.25;
end
if ~exist('weightMax','var')
    weightMax = 1;
end
if ~exist('winSize','var')
    winSize = -1;
end


if size(ts1, 2) == size(ts2, 2)
    dim = size(ts2, 2); %#ok<NASGU>
else
    error('Two time series should have the same dimension')
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
        tmp = TS_weightFun(abs(i - j), (m + n) / 2);
%         if ( tmp > 0.1 && tmp < 0.9)
%             fprintf('%d --- %d ---- %f ---- \n', i, j, tmp);
%         end
        
        tmpDist = TS_weightFun(abs(i - j), (m + n) / 2) * ...
            TS_euclideanDistance(ts1(i - 1, :), ts2(j - 1, :));
        minVal = min([D(i, j - 1), D(i - 1, j - 1), D(i - 1, j)]);
        D(i, j) = tmpDist + minVal;
        L(i, j) = min([L(i, j - 1), L(i - 1, j - 1), L(i - 1, j)]) + 1;

    end
end


wdtwDist = D(m + 1, n + 1);
wdtwLen = L(m + 1, n + 1);

function weightOut = TS_weightFun(difference, len)

weightOut = weightMax / (1 + exp(-g * (difference - (len / 2))));
% if ( weightOut > 0.4 && weightOut < 0.6)
%     fprintf('%f --- %f ---- %f ---- \n', difference, len, tmp);
% end

end
end

   