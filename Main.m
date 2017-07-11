clc
clear
close all 

global root
root = 'C:\Users\Amir\Desktop\github';

folder = TS_getAddress;


%% compute distance matrix for all datasets and savve to folder.distance

distanceList = 1:16; % the number associated to the distance functions
datasetList = 1:32; % the number that associated to datasets

TS_getDistanceForDataset(28, distanceList)




[all1NN, best1NN] = TS_do1NN(10, 5);