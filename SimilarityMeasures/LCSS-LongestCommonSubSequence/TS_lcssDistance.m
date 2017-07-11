function lcssDist = TS_lcssDistance(ts1, ts2, eps, winSize)

%     ************
%
%     Description
%     ----------
%     Compute the Longest-Common-Subsequence distance between time series ts1 and ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Options
%     ----------
%     opt eps       : double, a threshold for considering distance
%                             default value = Inf
%     opt winSize   : integer, temporal constraint on the warping window
%                              size. default value = -1
%
%     Returns
%     -------
%     lcssDist      : double, The Longest-Common-Subsequence distance between time series ts1 and ts2
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
if ~exist('eps','var')
    eps = inf;
end

if size(ts1, 2) == size(ts2, 2)
    dim = size(ts2, 2); %#ok<NASGU>
else
    error('Two time series dimension must be the same')
end

m = size(ts1, 1);
n = size(ts2, 1);

D = zeros(m + 1, n + 1);

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
        
        tmp = TS_euclideanDistance(ts1(i - 1, :), ts2(j - 1, :));
        if ( tmp < eps)
            D(i, j) = D(i - 1, j - 1) + 1;
        else
            D(i, j) = max([D(i, j - 1), D(i - 1, j)]);
        end
    end
end

lcssDist = 1 - D(m + 1, n + 1) / min([m , n]);
