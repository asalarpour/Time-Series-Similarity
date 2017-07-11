function  [all1NN, best1NN] = TS_do1NN(iterNum, foldNum)

%     ************
%
%     Description
%     ----------
%     return the 1-nn results with foldNum number of folds for cross
%     validation for iterNum number of iterations. for best variant of each
%     distance function and all variant of each distance function and all
%     datasets
% 
%     Parameters
%     ----------
%     param iterNum 	:  int, the number of replicants for 1-nn proces
%     param foldNum 	:  int, the number of folds for k-fold cross validation
%
%     Options
%     ----------
%
%     Returns
%     -------
%     all1NN            :  struct, the complete list of different distances with all variants, 1-nn 
%                               and k-fold cross validation with iterations for all dataset
%     best1NN           :  struct, the best result select based on 50% of train fold 
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

% set Data

allDatasetTS = cell(0);
folder = TS_getAddress();

for datasetNum = 1:32
    
    disp(datasetNum)
    [tsSet, tsSpec] = TS_selectDataset(datasetNum, folder);

    datasetFilesAdd = TS_getDistanceFilesAddress(datasetNum, folder);

    [allResult1NN, bestResult1NN] = TS_getResult1NN(datasetFilesAdd, tsSet, iterNum, foldNum);

    allDatasetTS{datasetNum}.folder = folder;
    allDatasetTS{datasetNum}.tsSpec = tsSpec;

    allDatasetTS{datasetNum}.allResult1NN = allResult1NN;
    allDatasetTS{datasetNum}.bestResult1NN = bestResult1NN;

end


all1NN = allDatasetTS{1}.allResult1NN;
best1NN = allDatasetTS{1}.bestResult1NN;

for i = 2: length(allDatasetTS)
    all1NN = [all1NN allDatasetTS{i}.allResult1NN]; %#ok<AGROW>
    best1NN = [best1NN allDatasetTS{i}.allResult1NN]; %#ok<AGROW>
end











        

