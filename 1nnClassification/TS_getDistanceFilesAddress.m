function datasetFilesAdd = TS_getDistanceFilesAddress(datasetNum, folder)

%     ************
%
%     Description
%     ----------
%     return the list of distance files addresses that saved in
%     folder.distance address and associated with datasetNum
% 
%     Parameters
%     ----------
%     param datasetNum  :  int, the number of MTS dataset between 1 to 32
%     param folder      :  struct, consist the needed address for dataset and associated distance files
%
%     Options
%     ----------
%
%     Returns
%     -------
%     datasetFilesAdd   :  struct, list of distance files address for specified datasest
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
%     December 2016 : Last revision: 03-Feb-2017
%     
%     ************

% extract distances files name
filesInFolder = dir(folder.distance);

datasetFilesAdd = struct('address', 0,'datasetNum', 0,  'distanceNum', 0, 'distanceName', 0,...
    'param1Idx', 0, 'param1Name', 0, 'param2Idx', 0, 'param2Name', 0);

count = 0;
for i = 1: size( filesInFolder, 1)
    tmpName = filesInFolder(i).name;
    tmpComp = ['Dataset_' num2str(datasetNum) '_'];
    if length(tmpName) > length(tmpComp)
        if strcmp(tmpComp, tmpName(1:length(tmpComp)))
            count = count + 1;
            datasetFilesAdd(count).address = [filesInFolder(i).folder '\' tmpName];
            uIdx = find(tmpName == '_');
            datasetFilesAdd(count).datasetNum = str2num(tmpName(uIdx(1) + 1: uIdx(2)-1)); %#ok<ST2NM>
            datasetFilesAdd(count).distanceNum = str2num(tmpName(uIdx(2) + 1: uIdx(3)-1)); %#ok<ST2NM>
            datasetFilesAdd(count).distanceName = tmpName(uIdx(3) + 1: uIdx(4)-1);

            
            if (length(uIdx) > 4)
                datasetFilesAdd(count).param1Idx = str2num(tmpName(uIdx(5) + 1: uIdx(6)-1)); %#ok<ST2NM>
                datasetFilesAdd(count).param1Name = tmpName(uIdx(4) + 1: uIdx(5)-1);
            end
            if (length(uIdx) > 6)
                datasetFilesAdd(count).param2Idx = str2num(tmpName(uIdx(7) + 1: uIdx(8)-1)); %#ok<ST2NM>
                datasetFilesAdd(count).param2Name = tmpName(uIdx(6) + 1: uIdx(7)-1);
            end
        end
    end
end

[~, sortIdx] = sort([datasetFilesAdd.distanceNum]);
datasetFilesAdd = datasetFilesAdd(sortIdx);