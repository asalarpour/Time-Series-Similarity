function erpDist = TS_erpDistance(ts1, ts2, gap, winSize)

%     ************
%
%     Description
%     ----------
%     Compute the Edit distance with Real Penalty between time series ts1 and ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Options
%     ----------
%     opt gap       :  1 x dim, a sample point as a refrence to claculate
%                               the penalty
%     opt winSize   : integer, temporal constraint on the warping window
%                              size. default value = -1
%
%     Returns
%     -------
%     erpDist       :  double, The Edit distance with Real Penalty between time series ts1 and ts2.
%
%     Other m-files required    : TS_euclideanDistance
%     Subfunctions              : none
%     MAT-files required        : none
%     
%     References
%     ----------
%     @inproceedings{chen2004marriage,
%       title={On the marriage of lp-norms and edit distance},
%       author={Chen, Lei and Ng, Raymond},
%       booktitle={Proceedings of the Thirtieth international conference on Very large data bases-Volume 30},
%       pages={792--803},
%       year={2004},
%       organization={VLDB Endowment}
%     }
%
%     Author
%     ----------
%     Amir Salarpour, Ph.D. Candidate, Artificial Intelligence
%     Bu-Ali Sina University, Hamedan, Iran, Dept. of Computer Engineering
%     email address : amir.salarpour@gmail.com  
%     Website       : http://www.salarpour.com
%     December 2016 : Last revision: 28-Jan-2017
%     
%     ************

if ~exist('winSize','var')
    winSize = -1;
end


if size(ts1, 2) == size(ts2, 2)
    dim = size(ts2, 2); 
else
    error('Two time series should have the same dimension')
end
if ~exist('gap','var')
    gap = zeros(1, dim);
end

m = size(ts1, 1);
n = size(ts2, 1);

D = zeros(m + 1, n + 1);

for i = 1: m
    D(i + 1, 1) = D(i, 1) + TS_euclideanDistance(gap, ts1(i , :));
end

for j = 1: n
    D(1, j + 1) = D(1, j) + TS_euclideanDistance(gap, ts2(j , :));
end

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
        tmpErp0 = D(i - 1, j) + TS_euclideanDistance(ts1(i - 1, :), gap);
        tmpErp1 = D(i, j - 1) + TS_euclideanDistance(gap, ts2(j - 1, :));
        tmpErp01 = D(i - 1, j - 1) + TS_euclideanDistance(ts1(i - 1, :), ts2(j - 1, :));
        D(i, j)  = min([tmpErp0, tmpErp1, tmpErp01]);
    end
end

erpDist = D( m + 1, n + 1);

end
