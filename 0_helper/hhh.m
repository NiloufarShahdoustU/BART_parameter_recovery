% this is for parameter recovery for BART TD model

% Set Up
close all;
sigBHFEEGcorr = struct();
unitIDs_Pts = {};
clusterTbl = [];
sigContacts_allPts = false(0,4);
moreImpulsive = false(0,1);
MNItrodeLocs_allPTs = [];
trodeLabels_allPts = {''};
anatomicalLocs_allPts = {''};
ptIDs_all = {''};
leftHemiContacts_allPts = [];
noWMpts = {};
[ptArray,bhvStruct,hazEEG] = BARTnumbers;
nPts = length(ptArray);
set(0,'defaultfigurerenderer','painters')
pts = 1:nPts;
ptMap = cbrewer('qual','Set3',nPts);

% Loading in TDdata for each patient
mainPath = '\\155.100.91.44\D\Data\preProcessed\BART_preprocessed\';
neuralData_All = {};

% Loading NeuralData All
for p = 1:length(ptArray)
    ptID = ptArray{p};
    dataPath = fullfile(mainPath, ptID, 'Data');
    fileInfo = dir(fullfile(dataPath, '*TDdata.mat'));

    if ~isempty(fileInfo)
        fileToLoad = fullfile(dataPath, fileInfo(1).name);
        learningrate_data = load(fileToLoad);
        neuralData_All{p} = learningrate_data;
    else
        warning('No TDdata.mat file found for %s', ptID);
        neuralData_All{p} = [];
    end
end

% Getting learning rates
for k = 1:nPts
    % Getting positive an negative learning rates
    rstd_positivePE(k) = neuralData_All{k}.TDdata.positiveBetaMaxIdx/100; % positivePE dividing by 100 because it needs to go from 18 to 0.18
    rstd_negativePE(k) = neuralData_All{k}.TDdata.negativeBetaMaxIdx/100; % negativePE
    % rstd metric: ratio of pos/neg
    RSTD_metric_neural(k) = (rstd_negativePE(k) - rstd_positivePE(k)) / (rstd_negativePE(k) + rstd_positivePE(k)); % calculating ratio.
end
