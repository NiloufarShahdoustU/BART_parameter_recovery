clear;
clc;
close all;
warning('off','all')





outputFolderName = 'param_recovery_4_param_recovery';

inputFolderName = '\\155.100.91.44\d\Data\Nill\BART_param_recovery\matlab_modeling\param_recovery_3_simulated_fields';

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName);
end


matFiles = dir(fullfile(inputFolderName, '*.mat'));
nPatients = length(matFiles);



ptID_all = cell(nPatients,1);
sim_alpha_plus_all = zeros(nPatients,1);
sim_alpha_minus_all = zeros(nPatients,1);
fit_alpha_plus_all = zeros(nPatients,1);
fit_alpha_minus_all = zeros(nPatients,1);


ptID_all = cell(nPatients,1);
for pt = 1:nPatients
% for pt = 1:1

    fileName = matFiles(pt).name;
    fprintf('processing pt: %s\n', fileName);

    matFilePath = fullfile(inputFolderName, fileName);
    load(matFilePath);
    nTrials = TDdataParamRecovery.nTrials;

    [~, ptID, ~] = fileparts(fileName);
    ptID = strrep(ptID, '_TDdataParamRecovery', '');


    result = TDdataParamRecovery.resultSimulated(1:nTrials);
    prevReward = TDdataParamRecovery.Reward(1:nTrials);

    outcomeType = nan(1,nTrials);
    outcomeType(strcmp(result,'banked')) = 1;
    outcomeType(strcmp(result,'popped')) = 2;

    % alpha grid
    a = 0.01:0.01:1;
    nAlpha = length(a);

    % store results for this patient
    inverseTemperatureRSTD = nan(nAlpha, nAlpha);
    rstdV = nan(nAlpha, nAlpha, nTrials);
    rstdRPE = nan(nAlpha, nAlpha, nTrials);
    points = TDdataParamRecovery.points;

    pointsPerTrial = TDdataParamRecovery.rewardSimulated;
    
    % now fitting risk-sensitive (asymmetric) models.

    Reward = zeros(1,nTrials);
    RewardPE = zeros(1,nTrials);
    Reward(1) = points(1);
    expectedReward = zeros(1,nTrials);
    a = 0.01:0.01:1;		% alphas;
    
    for ap = length(a):-1:1
        for an = length(a):-1:1
            for t = 2:nTrials
                % updating reward variable in each trial.
                if strcmp(result(t),'banked')
                    Reward(t) = pointsPerTrial(t);			% outcome on current trial.
                else
                    Reward(t) = 0;
                end
    
                % updating risk and reward PE
                RewardPE(t) = Reward(t) - expectedReward(t-1);
    
                % updating value in a piece wise fashion with different learning rates for postive and negative RPEs, respectively.
                if RewardPE(t)>0
                    expectedReward(t) = expectedReward(t-1) + a(ap)*RewardPE(t);% was t-1, this was wrong.
                elseif RewardPE(t)<0
                    expectedReward(t) = expectedReward(t-1) + a(an)*RewardPE(t);
                else
                    expectedReward(t) = expectedReward(t-1); % in the case when expectedReward == 0, expectedReward is expectedReward(t-1).
                end
            end
    
            whichLink = 'logit';
            % similar but using probit instead of logit. Logit gives the same
            % answer as above.
            B = glmfit(expectedReward',categorical(outcomeType),'binomial','link',whichLink);
    
            % inverse temperature parameters from glmfit
            inverseTemperatureRSTD(ap,an) = B(2);
    
            % TD vars
            rstdV(ap,an,:) = expectedReward;
            rstdRPE(ap,an,:) = RewardPE;
        end
    end
    
    % Getting gradient alphas (behavior)
    [~,positiveBetaMaxIdx] = max(max(inverseTemperatureRSTD,[],2),[],1);
    [~,negativeBetaMaxIdx] = max(max(inverseTemperatureRSTD,[],1),[],2);

    bestAlphaPos = a(positiveBetaMaxIdx);
    bestAlphaNeg = a(negativeBetaMaxIdx);

    sim_alpha_plus = a(positiveBetaMaxIdx);
    sim_alpha_minus = a(negativeBetaMaxIdx);

    fit_alpha_plus = TDdataParamRecovery.bestAlphaPos;
    fit_alpha_minus = TDdataParamRecovery.bestAlphaNeg;

    ptID_all{pt} = ptID;
    sim_alpha_plus_all(pt) = sim_alpha_plus;
    sim_alpha_minus_all(pt) = sim_alpha_minus;
    fit_alpha_plus_all(pt) = fit_alpha_plus;
    fit_alpha_minus_all(pt) = fit_alpha_minus;
end

T = table(ptID_all, sim_alpha_plus_all, sim_alpha_minus_all, ...
          fit_alpha_plus_all, fit_alpha_minus_all, ...
          'VariableNames', {'ptID','sim_alpha_plus','sim_alpha_minus','fit_alpha_plus','fit_alpha_minus'});

writetable(T, fullfile(outputFolderName, 'alpha_comparison.csv'));




