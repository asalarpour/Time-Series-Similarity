function [allResult1NN, bestResult1NN] = TS_getResult1NN(datasetFilesAdd, tsSet, iterNum, foldNum)

%     ************
%
%     Description
%     ----------
%     datasetFilesAdd consist the address of all calculated distances for
%     specific dataset (tsSet) the 1-NN with the foldNum number of folds is
%     done and iterated for iterNum times 
% 
%     Parameters
%     ----------
%     param datasetFilesAdd :  struct, the list of distance file addresses for specified dataset
%     param tsSet           :  struct, the specified dataset
%     param iterNum         :  int, the number of replicants for 1-nn proces
%     param foldNum         :  int, the number of folds for k-fold cross validation
%
%     Options
%     ----------
%
%     Returns
%     -------
%     allResult1NN          :  struct, the complete list of different distances with all variants, 1-nn 
%                               and k-fold cross validation with iterations for specific dataset
%     bestResult1NN         :  struct, the best result select based on 50% of train fold 
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



obsNum = length(tsSet);
obs = (1:obsNum)';
label = [tsSet(:).label]';

reverseStr = '';

idxFold = [];
for i = 1:iterNum
    idxFold = [idxFold crossvalind('Kfold', label, foldNum)]; %#ok<AGROW>
end

allResult1NN = struct('datasetNum', [], 'distNum', [], 'distName', [],...
    'param1Idx', [], 'param1Name', [], 'param2Idx', [], 'param2Name', [],...
    'time',[], 'num', [], 'accLvo', [], 'accTst', [] );

if all(datasetFilesAdd(1).address == 0)
    bestResult1NN = allResult1NN;
    return
end

for i = 1:length(datasetFilesAdd)
    
    % Display the progress
    percentDone = 100 * i / length(datasetFilesAdd) ;
    msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    allResult1NN(i).datasetNum = datasetFilesAdd(i).datasetNum;
    allResult1NN(i).distNum = datasetFilesAdd(i).distanceNum;
    allResult1NN(i).distName = datasetFilesAdd(i).distanceName;
    allResult1NN(i).param1Idx = datasetFilesAdd(i).param1Idx;
    allResult1NN(i).param1Name = datasetFilesAdd(i).param1Name;
    allResult1NN(i).param2Idx = datasetFilesAdd(i).param2Idx;
    allResult1NN(i).param2Name = datasetFilesAdd(i).param2Name;    
    
    allResult1NN(i).num = i;
    
    load(datasetFilesAdd(i).address);
    
    allResult1NN(i).time = distCell{1}.totalTime;
    distSrc = squareform(distCell{1}.distMat);
    
    if (allResult1NN(i).distNum == 6)
        distSrc = squareform( 1 - (distCell{1}.distMat ./ max(distCell{1}.distMat)));
    end
    
    accLvo = [];
    accTst = [];
    for iter = 1: iterNum
        idxTmp = idxFold(:, iter);
        
        lvoTmp = [];
        tstTmp = [];
        
        for fold = 1: foldNum
            
            trnObs = obs(idxTmp == fold);
            trnLbl = label(idxTmp == fold);
            tstObs = obs(idxTmp ~= fold);
            tstLbl = label(idxTmp ~= fold);
            
            distTmp = distSrc(trnObs, trnObs);
            distTmp = distTmp + eye(size(distTmp)) * max(max(distTmp));
            [~, minIdxTmp] = min(distTmp);
            lvoAcc = sum(trnLbl(minIdxTmp) == trnLbl) / length(trnLbl);
            
            distTmp = distSrc(trnObs, tstObs);
            if size(distTmp , 1) == 1
                tstAcc = sum(trnLbl == tstLbl) / length(tstLbl);
            else
                [~, minIdxTmp] = min(distTmp);
                tstAcc = sum(trnLbl(minIdxTmp) == tstLbl) / length(tstLbl);
            end          
            
            lvoTmp = [lvoTmp lvoAcc]; %#ok<AGROW>
            tstTmp = [tstTmp tstAcc]; %#ok<AGROW>
            
            
        end
        accLvo = [accLvo; lvoTmp]; %#ok<AGROW>
        accTst = [accTst; tstTmp]; %#ok<AGROW>
    end
    allResult1NN(i).accLvo = accLvo;
    allResult1NN(i).accTst = accTst;
end
    

tmpDist = [allResult1NN(:).distNum];
tmpUnq = unique(tmpDist);

for i = 1: length(tmpUnq)
    tmpIdx = find([allResult1NN(:).distNum] == tmpUnq(i));
    if (length(tmpIdx) == 1)
        
        tmpStruct = allResult1NN(tmpIdx); 
        tmpStruct.accMean = mean(mean(tmpStruct.accTst)); 
        tmpStruct.allLvo = []; 
        tmpStruct.allTst = []; 
        tmpStruct.bstLvoList = []; 
        
        bestResult1NN(i) = tmpStruct; 
    else
        allLvo = [];
        allTst= [];
        for j = 1: length(tmpIdx)
            allLvo = [allLvo mean(allResult1NN(tmpIdx(j)).accLvo, 2)]; %#ok<AGROW>
            allTst = [allTst mean(allResult1NN(tmpIdx(j)).accTst, 2)]; %#ok<AGROW>
        end
        
        [~, bstLvoList] = max(allLvo');
        
        tmpBst = mode(bstLvoList);
        
        tmpStruct = allResult1NN(tmpIdx(tmpBst)); 
        tmpStruct.accMean = mean(mean(tmpStruct.accTst)); 
        tmpStruct.allLvo = allLvo; 
        tmpStruct.allTst = allTst; 
        tmpStruct.bstLvoList = bstLvoList; 
        
        bestResult1NN(i) = tmpStruct; 
    end
end

