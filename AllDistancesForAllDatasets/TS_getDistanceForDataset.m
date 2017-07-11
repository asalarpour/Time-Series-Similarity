function TS_getDistanceForDataset(datasetList, distanceList)

% dist 1  -----> DTW Full
% dist 2  -----> DTW Constrained
% dist 3  -----> LCSS EPSILON
% dist 4  -----> EDR EPSILON
% dist 5  -----> ERP GAP
% dist 6  -----> SWALE EPSILON GAP
% dist 7  -----> TWED LAMBDA NU
% dist 8  -----> MSM COST
% dist 9  -----> WDTW G
% dist 10 -----> Hausdorff
% dist 11 -----> frechet
% dist 12 -----> discreteFrechet
% dist 13 -----> sspd
% dist 14 -----> TQuEST THRESHOLD
% dist 15 -----> CIDTW
% dist 16 -----> DDTW ALPHA

for datasetNum = datasetList
    for distNum = distanceList
        if  distNum == 14 || distNum == 11 % distance 14 and 11 are eliminated for MTS datasets
            continue
        end

        TS_computeDistanceforDatasetAllInMex(datasetNum, distNum, 5)
    end
end

