clear;
clc;
close all;


inputFolderName  = '\\155.100.91.44\d\Data\Nill\BART_param_recovery\param_recovery_1_data\';

%%

matFiles = dir(fullfile(inputFolderName, '*.mat'));
nPatients = length(matFiles);

% for pt = 1:nPatients
for pt = 3:3   
    fileName = matFiles(pt).name;
    [~, ptID, ~] = fileparts(fileName);
    ptID = erase(ptID, '_TDdataParamRecovery');

    fprintf('processing pt %d/%d: %s\n', pt, nPatients, ptID);

    matFile = fullfile(inputFolderName, fileName);
    load(matFile);
end

%% taking a look at the points - Reward data

results = TDdataParamRecovery.result;
isControl = TDdataParamRecovery.is_control;

%the control trials don't pop
pointsMinusRewardPop  = TDdataParamRecovery.pointsMinusReward(strcmp(results,'popped'));

pointsMinusRewardBankNotControl = TDdataParamRecovery.pointsMinusReward(strcmp(results,'banked')& isControl == 0);
pointsMinusRewardBankControl = TDdataParamRecovery.pointsMinusReward(strcmp(results,'banked')& isControl == 1);

%% taking a look at the trial type for all patients:
maxTrialNum = 300;
trialTypeAcrossPts = nan(maxTrialNum, nPatients);
ptIDs = cell(nPatients,1);

for pt = 1:nPatients

    fileName = matFiles(pt).name;
    [~, ptID, ~] = fileparts(fileName);
    ptID = erase(ptID, '_TDdataParamRecovery');

    ptIDs{pt} = ptID;

    fprintf('processing pt %d/%d: %s\n', pt, nPatients, ptID);

    matFile = fullfile(inputFolderName, fileName);
    load(matFile);

    trialType = TDdataParamRecovery.trial_type;
    nTrials = length(trialType);

    trialTypeAcrossPts(1:nTrials, pt) = trialType;

end


% Make columns equal length

minTrials = min(sum(~isnan(trialTypeAcrossPts)));
trialTypeAcrossPtsSamesize = trialTypeAcrossPts(1:minTrials, :);


% Sort columns so identical vectors are next to each other

[~, order] = sortrows(trialTypeAcrossPtsSamesize');
trialTypeGrouped = trialTypeAcrossPtsSamesize(:, order);


% Find unique trial-type vectors

[uniqueCols, ~, ic] = unique(trialTypeAcrossPtsSamesize','rows');

uniquecounts = accumarray(ic,1);


% Print ptIDs belonging to each unique vector

for i = 1:size(uniqueCols,1)

    pts = ptIDs(ic == i);

    fprintf('\nTrial type vector %d (count = %d):\n', i, uniquecounts(i));
    disp(pts)

end