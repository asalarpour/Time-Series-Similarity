function [freeSpace1, freeSpace2] = TS_freeSpace(ts1, ts2, eps)

%     ************
%
%     Description
%     ----------
%     Compute all the free space on the boundary of cells in the diagram for polygonal chains ts1 and ts1 and the given eps
%     freeSpace1{i,j} is the free space of segment [ts1(i, :),ts1(i + 1, :)] from point ts2(j, :)
%     freeSpace1{i,j} is the free space of segment [ts2(j, :),ts2(j + 1, :)] from point ts1(i, :)
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%     param eps     :  double, reachability distance
%
%     Returns
%     -------
%     freeSpace1    : cell, free spaces of segments of ts1 from points of ts2
%     freeSpace2    : cell, free spaces of segments of ts2 from points of ts1
%
%     Other m-files required    : TS_freeLine
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

freeSpace1 = cell(1);
for j = 1: n
    for i = 1: m - 1
        freeSpace1{i, j} = TS_freeLine(ts2(j, :), eps, ts1(i: i + 1, :));
    end
end


freeSpace2 = cell(1);
for j = 1:  n - 1
    for i = 1: m
        freeSpace2{i, j} = TS_freeLine(ts1(i, :), eps, ts2(j: j + 1, :));
    end
end
end
