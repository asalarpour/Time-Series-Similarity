function datasetsInfo = TS_getDatasetInfo()

%     ************
%
%     Description
%     ----------
%     return the time series dataset information
% 
%     Parameters
%     ----------
%
%     Options
%     ----------
%
%     Returns
%     -------
%     datasetsInfo       :   struct, return the needed information.
%
%     Other m-files required    : TS_getAddress, TS_euclideanDistance
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

datasetsInfo = struct('name', 'a', 'tsCount', 0, 'classCount', 0, 'complexityAvg', 0, ...
     'pointCountSTD', 0, 'tsLenSTD', 0);

reverseStr = ''; 
 
for k = 1: 32
    
    % Display the progress
    percentDone = 100 * k / 32 ;
    msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    [tsSet, tsSpec] = TS_selectDataset(k, folder);

    datasetsInfo(k).tsCount = length(tsSet);
    datasetsInfo(k).name = tsSpec.name;
    
    datasetsInfo(k).classCount = length(unique([tsSet.label]));

    tsInfo = struct('pointNum', 0, 'tsLen', 0, 'lenDirect', 0, 'complexity', 0);

    for i = 1: datasetsInfo(k).tsCount
        tmp = tsSet(i).ts;
        tmpLen = 0;
        for j = 1: size(tmp, 1) - 1
            tmpLen = tmpLen + TS_euclideanDistance(tmp(j, :), tmp(j + 1, :));
        end
        tmpFirst2Final = TS_euclideanDistance(tmp(1, :), tmp(end, :));
        tsInfo(i).pointNum = size(tmp, 1);
        tsInfo(i).tsLen = tmpLen;
        tsInfo(i).lenDirect = tmpFirst2Final;
        tsInfo(i).complexity = tmpFirst2Final / tmpLen;
    end

    datasetsInfo(k).dim = size(tsSet(1).ts, 2);
    datasetsInfo(k).complexityAvg = mean([tsInfo.complexity]);
    datasetsInfo(k).meanPoints = mean([tsInfo.pointNum]);
    
    datasetsInfo(k).pointCountSTD = std([tsInfo.pointNum]);
    datasetsInfo(k).maxPoints = max([tsInfo.pointNum]);
    datasetsInfo(k).minPoints = min([tsInfo.pointNum]);
    
    datasetsInfo(k).tsLenSTD = std([tsInfo.tsLen]);
    datasetsInfo(k).tsLenAVG = mean([tsInfo.tsLen]);
    datasetsInfo(k).tsLenMAX = max([tsInfo.tsLen]);
    datasetsInfo(k).tsLenMIN = min([tsInfo.tsLen]);
    

end

