function ddtwDist = TS_deriviativeDtwDistanceMex(ts1, ts2, alpha)
    
%     ************
%
%     Description
%     ----------
%     Compute the Deriviative Dynamic-Time Warping distance between time series ts1 and ts2.
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
%     ddtwDist       : double, The Deriviative Dynamic-Time Warping distance between time series ts1 and ts2
%
%     Other m-files required    : TS_dtwDistanceMex
%     Subfunctions              : none
%     MAT-files required        : none
%     
%     References
%     ----------
%     @article{gorecki2013using,
%       title={Using derivatives in time series classification},
%       author={G{\'o}recki, Tomasz and {\L}uczak, Maciej},
%       journal={Data Mining and Knowledge Discovery},
%       pages={1--22},
%       year={2013},
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

if ~exist('alpha','var')
    alpha = 1;
end


a = cos(alpha);
b = sin(alpha);

dist1 = TS_dtwDistanceMex(ts1, ts2);
dist2 = TS_dtwDistanceMex(diff(ts1, 1, 1), diff(ts2, 1, 1));

ddtwDist = a * dist1 + b * dist2;

end

   