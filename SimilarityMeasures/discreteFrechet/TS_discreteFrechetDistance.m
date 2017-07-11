function disFrechDist = TS_discreteFrechetDistance( ts1, ts2)

%     ************
%
%     Description
%     ----------
%     Compute the discret frechet distance between time series ts1 and ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Returns
%     -------
%     disFrechDist  :  double, the discret frechet distance between time series ts1 and ts2
%
%     Other m-files required    : TS_euclideanDistance
%     Subfunctions              : TS_func
%     MAT-files required        : none
%     
%     References
%     ----------
%     @techreport{eiter1994computing,
%       title={Computing discrete Fr{\'e}chet distance},
%       author={Eiter, Thomas and Mannila, Heikki},
%       year={1994},
%       institution={Citeseer}
%     }
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
ca = -ones(m, n);
disFrechDist =  TS_func(m, n, ts1, ts2);
    


% =======================================================================
function output = TS_func(i, j, ts1, ts2)

if (ca(i, j) > -1)
    output = ca(i, j);
    return
    
elseif (i == 1 && j == 1)
    ca(i, j) = TS_euclideanDistance(ts1(1, :), ts2(1, :));
    
elseif (i > 1 && j == 1)
    ca(i, j) = max([TS_func(i - 1, 1, ts1, ts2),...
        TS_euclideanDistance(ts1(i, :),ts2(1, :))]);
    
elseif (i == 1 && j > 0)
    ca(i,j) = max([TS_func(1, j - 1, ts1, ts2),...
        TS_euclideanDistance(ts1(1, :), ts2(j, :))]);
    
elseif (i > 1 && j > 1)
%     disp([num2str(i), '------', num2str(j)])

    ca(i,j) = max([min([TS_func(i - 1, j, ts1, ts2),...
        TS_func(i - 1, j - 1, ts1, ts2),...
        TS_func(i, j - 1, ts1, ts2)]),...
        TS_euclideanDistance(ts1(i, :),ts2(j, :))]);
else
    ca(i, j) = inf;
end
    output = ca(i,j);
end
end
    
    
    
    
    
    
    
    
    
    