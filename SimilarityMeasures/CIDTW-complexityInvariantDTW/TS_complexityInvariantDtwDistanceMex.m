function cidtwDist = TS_complexityInvariantDtwDistanceMex(ts1, ts2)
    
%     ************
%
%     Description
%     ----------
%     Compute the Complexity Invariant Dynamic-Time Warping distance between time series ts1 and ts2.
% 
%     Parameters
%     ----------
%     param ts1   :  m x dim, time series 1 matrix with the length of m
%     param ts2   :  n x dim, time series 2 matrix with the length of n
%
%     Options
%     ----------
%     opt alpha     :  double, the test interval is between 1 and pi / 2 with 0.01 step
%
%     Returns
%     -------
%     cidtwDist     :  double, The Complexity Invariant Dynamic-Time Warping distance between time series ts1 and ts2
%
%     Other m-files required    : TS_dtwDistance
%     Subfunctions              : none
%     MAT-files required        : none
%     
%     References
%     ----------
%     @article{batista2014cid,
%       title={CID: an efficient complexity-invariant distance for time series},
%       author={Batista, Gustavo EAPA and Keogh, Eamonn J and Tataw, Oben Moses and De Souza, Vinicius MA},
%       journal={Data Mining and Knowledge Discovery},
%       volume={28},
%       number={3},
%       pages={634--669},
%       year={2014},
%       publisher={Springer}
%     }
%
%     Author
%     ----------
%     Amir Salarpour, Ph.D. Candidate, Artificial Intelligence
%     Bu-Ali Sina University, Hamedan, Iran, Dept. of Computer Engineering
%     email address : amir.salarpour@gmail.com  
%     Website       : http://www.salarpour.com
%     December 2016 : Last revision: 03-Feb-2017
%     
%     ************


dist1 = TS_dtwDistanceMex(ts1, ts2);

CE1 = sqrt(sum(sum(diff(ts1, 1,1) .^ 2)));
CE2 = sqrt(sum(sum(diff(ts2, 1,1) .^ 2)));


cidtwDist = dist1 * (max([CE1, CE2]) / min([CE1, CE2]));

end

   