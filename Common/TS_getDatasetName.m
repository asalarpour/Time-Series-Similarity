function dataName = TS_getDatasetName(DatasetNum)

%     ************
%
%     Description
%     ----------
%     return the time series dataset name
% 
%     Parameters
%     ----------
%
%     Options
%     ----------
%
%     Returns
%     -------
%     datasetsInfo       :   struct, return the dataset name.
%
%     Other m-files required    : TS_getAddress
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

folder = TS_getAddress();


dataName = cell(0);
for i = 1:DatasetNum
    [~, tsSpec] = TS_selectDataset(i, folder);
    dataName{i} = tsSpec.name;

end