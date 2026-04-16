clear;
clc;
close all;
warning('off','all')


%%

inputFolderName = '\\155.100.91.44\d\Data\Nill\BART_param_recovery\context_modeling\param_recovery_1_modeling\';



matFiles = dir(fullfile(inputFolderName, '*.mat'));
nPatients = length(matFiles);

allAlphaPos = zeros(1,nPatients);
allAlphaNeg = zeros(1,nPatients);

for pt = 1:nPatients
% for pt = 1:1

    fileName = matFiles(pt).name;
    fprintf('processing pt: %s\n', fileName);

    matFilePath = fullfile(inputFolderName, fileName);
    load(matFilePath);
    nTrials = TDdataParamRecovery.nTrials;

    bestAlphaPos = TDdataParamRecovery.bestAlphaPos;
    bestAlphaNeg = TDdataParamRecovery.bestAlphaNeg;

    allAlphaPos(pt) = bestAlphaPos;
    allAlphaNeg(pt) = bestAlphaNeg;

end

%%

figure;
scatter(allAlphaPos, allAlphaNeg, 20, 'filled');
xlabel('\alpha_+');
ylabel('\alpha_-');
