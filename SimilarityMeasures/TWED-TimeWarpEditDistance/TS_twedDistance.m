function twedDist = TS_twedDistance(ts1, ts2, lambda, nu, winSize)

%     ************
%
%     Description
%     ----------
%     Compute the Time Warp Edit distance between time series ts1 and ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Options
%     ----------
%     opt lambda    : double, penalty, punishment for distances at deletions
%     opt nu        : double, stiffness, determines the elasticity Nu > = 0 required for distance measurement.
%     opt winSize   : integer, temporal constraint on the warping window
%                              size. default value = -1
%
%     Returns
%     -------
%     twedDist      : double, The Time Warp Edit distance between time series ts1 and ts2 
%
%     Other m-files required    : TS_euclideanDistance
%     Subfunctions              : none
%     MAT-files required        : none
%     
%     References
%     ----------
%     http://people.irisa.fr/Pierre-Francois.Marteau/
%     @article{marteau2009time,
%       title={Time warp edit distance with stiffness adjustment for time series matching},
%       author={Marteau, Pierre-Fran{\c{c}}ois},
%       journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
%       volume={31},
%       number={2},
%       pages={306--318},
%       year={2009},
%       publisher={IEEE}
%     }
%
%     Author
%     ----------
%     Amir Salarpour, Ph.D. Candidate, Artificial Intelligence
%     Bu-Ali Sina University, Hamedan, Iran, Dept. of Computer Engineering
%     email address : amir.salarpour@gmail.com  
%     Website       : http://www.salarpour.com
%     December 2016 : Last revision: 26-Jan-2017
%     
%     ************

if ~exist('winSize','var')
    winSize = -1;
end
if ~exist('lambda','var')
    lambda = 1;
end
if ~exist('nu','var')
    nu = 0;
end

if size(ts1, 2) == size(ts2, 2)
    dim = size(ts2, 2); %#ok<NASGU>
else
    error('Two time series should have the same dimension')
end

if nu < 0
    warning('nu is negative')
    return
end

m = size(ts1, 1);
n = size(ts2, 1);

A = [zeros(1, size(ts1, 2)); ts1];
B = [zeros(1, size(ts2, 2)); ts2];

% This algorithm uses dynamic programming 
D = zeros(m + 1, n + 1);

% Initialize the D matrix, the first row and column set to infinitely
D(1, :) = inf;
D(:, 1) = inf;
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
        
        tmpUpdate1 = D(i - 1, j) + ...
            TS_euclideanDistance(A(i - 1, :), A(i, :)) + nu + lambda;
        
        tmpUpdate2 = D(i , j - 1) + ...
            TS_euclideanDistance(B(j - 1, :), B(j, :)) + nu + lambda;
        
        tmpUpdate3 = D(i - 1, j - 1) + ...
            TS_euclideanDistance(A(i, :), B(j, :)) + ...
            TS_euclideanDistance(A(i - 1, :), B(j - 1, :)) + (nu * 2 * abs(i - j)) ;
        
        D(i, j) = min([tmpUpdate1, tmpUpdate2, tmpUpdate3]);

    end
end

twedDist = D(m + 1, n + 1);

end