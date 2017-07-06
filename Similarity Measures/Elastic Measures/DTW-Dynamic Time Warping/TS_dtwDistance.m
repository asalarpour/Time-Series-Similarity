function [dtwDist, dtwLen] = TS_dtwDistance(traj1, traj2, winSize)
    
%     ************
%
%     Description
%     ----------
%     Compute the Dynamic-Time Warping distance between trajectory traj1 and traj2.
% 
%     Parameters
%     ----------
%     param traj1   :   m x dim, trajectory 1 matrix with the length of m
%     param traj2   :   n x dim, trajectory 2 matrix with the length of n
%
%     Options
%     ----------
%     opt winSize   :   integer,   temporal constraint on the warping window
%                              size. default value = -1
%
%     Returns
%     -------
%     dtwDist       :   double, The Dynamic-Time Warping distance between trajectory traj1 and traj2
%     dtwLen        :   integer, Show the length of warping path between trajectories
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

if size(traj1, 2) == size(traj2, 2)
    dim = size(traj2, 2); %#ok<NASGU>
else
    error('Two trajectories dimension must be the same')
end

m = size(traj1, 1);
n = size(traj2, 1);

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
        
        tmpDist = TS_euclideanDistance(traj1(i - 1, :), traj2(j - 1, :));
        minVal = min([D(i, j - 1), D(i - 1, j - 1), D(i - 1, j)]);
        D(i, j) = tmpDist + minVal;
        L(i, j) = min([L(i, j - 1), L(i - 1, j - 1), L(i - 1, j)]) + 1;

    end
end



dtwDist = D(m + 1, n + 1);
dtwLen = L(m + 1, n + 1);

end


   