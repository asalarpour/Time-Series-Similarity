function folder = TS_getAddress()

%     ************
%
%     Description
%     ----------
%     return the needed address for time series analysis
% 
%     Parameters
%     ----------
%     param root   :   the address of root folder
%
%     Options
%     ----------
%
%     Returns
%     -------
%     folder       :   struct, return the needed addresses.
%
%     Other m-files required    : none
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

global root

folder.dataset = [root '\TimeSeriesAnalysis\Data\MTS-Datasets\'];
folder.distance = [root '\TimeSeriesAnalysis\Data\MTS-Distance\'];
folder.figure = [root '\TimeSeriesAnalysis\Data\Figures\'];

addpath(genpath([root '\TimeSeriesAnalysis\']));
