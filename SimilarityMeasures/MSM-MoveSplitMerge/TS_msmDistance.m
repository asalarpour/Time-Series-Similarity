function msmDist = TS_msmDistance(ts1, ts2, cost, winSize)

%     ************
%
%     Description
%     ----------
%     Compute the Move Split Merge between time series ts1 and ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Options
%     ----------
%     opt cost      : double, cost of Split/Merge operation.
%     opt winSize   : integer, temporal constraint on the warping window
%                              size. default value = -1
%
%     Returns
%     -------
%     msmDist      : double, The Move Split Merge distance between time series ts1 and ts2 
%
%     Other m-files required    : TS_euclideanDistance
%     Subfunctions              : TS_msmCost
%     MAT-files required        : none
%     
%     References
%     ----------
%     @article{stefan2013move,
%       title={The move-split-merge metric for time series},
%       author={Stefan, Alexandra and Athitsos, Vassilis and Das, Gautam},
%       journal={IEEE transactions on Knowledge and Data Engineering},
%       volume={25},
%       number={6},
%       pages={1425--1438},
%       year={2013},
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
if ~exist('cost','var')
    cost = 0.1;
end

if size(ts1, 2) == size(ts2, 2)
    dim = size(ts2, 2); %#ok<NASGU>
else
    error('Two time series should have the same dimension')
end

m = size(ts1, 1);
n = size(ts2, 1);

D = zeros(m, n);

D(1, 1) = TS_euclideanDistance(ts1(1, :), ts2(1, :));

if winSize ~= -1
    winSize = max([winSize, abs(n - m)]);
end

for i = 2: m
    D(i, 1) = D(i - 1, 1) + TS_msmCost(ts1(i, :), ts1(i - 1, :), ts2(1, :), cost);
end

for j = 2: n
    D(1, j) = D(1, j - 1) + TS_msmCost(ts2(j, :), ts1(1, :), ts2(j - 1, :), cost);
end

for i = 2: m
    if winSize == -1
        jS = 2;
        jF = n;
    else
        jS = max([2, i - winSize]);
        jF = min ([n, i + winSize]);
    end    
    for j = jS:jF
        tmpDist1 = D(i - 1, j - 1) + TS_euclideanDistance( ts1(i, :), ts2(j, :));
        tmpDist2 = D(i - 1, j) + TS_msmCost(ts1(i, :), ts1(i - 1, :), ts2(j, :), cost);
        tmpDist3 = D(i, j - 1) + TS_msmCost(ts2(j, :), ts1(i, :), ts2(j - 1, :), cost);
        D(i, j) = min([tmpDist1, tmpDist2, tmpDist3]);
    end
end

msmDist = D(m, n);

end

function dist = TS_msmCost( newPoint, point1, point2, cost)

ax = newPoint(1);
ay = newPoint(2);
bx = point1(1);
by = point1(2);
cx = point2(1);
cy = point2(2);

if ( ( (bx <= ax) && (by <= ay) && (ax <= cx) &&(ay <= cy)) || ...
     ( (cx <= ax) && (cy <= ay) && (ax <= bx) && (ay <= by)) )
    dist = cost;
else
    dist = cost + min ( TS_euclideanDistance(newPoint, point1),...
        TS_euclideanDistance(newPoint, point2) );
end

end
