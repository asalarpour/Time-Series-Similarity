function [tsSet, tsSpec] = TS_selectDataset(datasetNum, folder)

%     ************
%
%     Description
%     ----------
%     return the specified time series dataset
% 
%     Parameters
%     ----------
%
%     Options
%     ----------
%	  datasetNum 		:	integer, the number of requested dataset
%	  folder 			:	struct, consist the related address 
%
%     Returns
%     -------
%     tsSet       		:   struct, the sample and label information
%	  tsSpec			:	struct, the spec of requested dataset
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

name = [folder.dataset, 'dataset_', num2str(datasetNum), '.mat'];
load(name)

tsSet = data;
tsSpec = spec;