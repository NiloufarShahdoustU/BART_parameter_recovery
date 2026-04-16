
clear;
clc;
close all;
%% loading gradient data

outputFolderName = '\\155.100.91.44\d\Data\Nill\BART_param_recovery\context_modeling\param_recovery_0_data\';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end

inputFolder = '\\155.100.91.44\d\Data\preProcessed\BART_preprocessed';

allItems = dir(inputFolder);
allItems = allItems([allItems.isdir]);
allItems = allItems(~ismember({allItems.name}, {'.', '..'}));

%  only folder names that start with '20'
folderMask = startsWith({allItems.name}, '20');
targetFolders = allItems(folderMask);
nPatients = length(targetFolders);

for pt = 1:nPatients
% for pt = 2:2

    ptID = targetFolders(pt).name;
    fprintf('processing pt: %s\n', ptID);

    dataFolder = fullfile(inputFolder, ptID, 'Data');
    matFile = dir(fullfile(dataFolder, '*TDdataGradients.mat'));

    if isempty(matFile)
        fprintf('No TDdataGradients.mat file found in: %s\n', dataFolder);
        continue;
    end

    matPath = fullfile(dataFolder, matFile(1).name);
    TDdata = load(matPath);

    % gathering required data:
    TDdataParamRecovery.a = TDdata.a;
    TDdataParamRecovery.nTrials = TDdata.nTrials;
    scoreVec = [TDdata.dataBHV.score];
    TDdataParamRecovery.scoreVec = scoreVec;
    pointsVec = [TDdata.dataBHV.points];
    Reward = [pointsVec(1), diff(scoreVec)];
    TDdataParamRecovery.Reward = Reward;
    TDdataParamRecovery.result = {TDdata.dataBHV.result};
    pointsMinusReward = pointsVec - Reward;
    TDdataParamRecovery.pointsMinusReward = pointsMinusReward;
    TDdataParamRecovery.inflate_time = [TDdata.dataBHV.inflate_time];
    TDdataParamRecovery.points = [TDdata.dataBHV.points];
    TDdataParamRecovery.is_control = [TDdata.dataBHV.is_control];
    TDdataParamRecovery.trial_type = [TDdata.dataBHV.trial_type];
    TDdataParamRecovery.RPE = TDdata.TDdataGradients.rstdRPE;
    TDdataParamRecovery.expectedReward = TDdata.TDdataGradients.rstdV;
    TDdataParamRecovery.bestAlphaPos = TDdata.TDdataGradients.neuralFit.bestAlphaPositive;
    TDdataParamRecovery.bestAlphaNeg = TDdata.TDdataGradients.neuralFit.bestAlphaNegative;

    save(fullfile(outputFolderName, [ptID '_TDdataParamRecovery.mat']), 'TDdataParamRecovery');

end


