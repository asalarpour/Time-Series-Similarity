function swaleDist = TS_swaleDistance(ts1, ts2, eps, gapC, reward, winSize)

%     ************
%
%     Description
%     ----------
%     Compute the Seuence Weighted Alignment distance between time series ts1 and ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Options
%     ----------
%     opt eps       : double, a threshold for considering distance
%     opt gapC      : double, gap cost
%     opt reward    : double, match reward
%     opt winSize   : integer, temporal constraint on the warping window
%                              size. default value = -1
%
%     Returns
%     -------
%     swaleDist       : double, The Seuence Weighted Alignment distance between time series ts1 and ts2
%
%     Other m-files required    : TS_euclideanDistance
%     Subfunctions              : none
%     MAT-files required        : none
%     
%     References
%     ----------
%     @inproceedings{morse2007efficient,
%       title={An efficient and accurate method for evaluating time series similarity},
%       author={Morse, Michael D and Patel, Jignesh M},
%       booktitle={Proceedings of the 2007 ACM SIGMOD international conference on Management of data},
%       pages={569--580},
%       year={2007},
%       organization={ACM}
%     }   
%
%     Author
%     ----------
%     Amir Salarpour, Ph.D. Candidate, Artificial Intelligence
%     Bu-Ali Sina University, Hamedan, Iran, Dept. of Computer Engineering
%     email address : amir.salarpour@gmail.com  
%     Website       : http://www.salarpour.com
%     December 2016 : Last revision: 31-Jan-2017
%     
%     ************

if ~exist('winSize','var')
    winSize = -1;
end
if ~exist('eps','var')
    eps = inf;
end
if ~exist('gapC','var')
    gapC = 0;
end
if ~exist('reward','var')
    reward = 0;
end


if size(ts1, 2) == size(ts2, 2)
    dim = size(ts2, 2);  %#ok<NASGU>
else
    error('Two time series should have the same dimension')
end


m = size(ts1, 1);
n = size(ts2, 1);

D = zeros(m + 1, n + 1);

D(1, :) = gapC * (0: n);
D(:, 1) = gapC * (0: m);

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
            D(i, j) = D(i - 1, j - 1) + reward;
        else
            D(i, j) = max([gapC + D(i, j - 1), gapC + D( i - 1, j)]);
        end
    end
end

swaleDist = D(m + 1, n + 1);

end